def test_embeddings(
    default_vllm_config, dist_init, num_loras, device, vocab_size, stage
) -> None:
    # For multi-GPU testing of Triton kernel, we must explicitly set the CUDA
    # device, see: https://github.com/triton-lang/triton/issues/2925
    # Same below.
    if current_platform.is_cuda_alike() or current_platform.is_xpu():
        torch.accelerator.set_device_index(device)

    torch.set_default_device(device)
    max_loras = 8
    lora_config = LoRAConfig(
        max_loras=max_loras, max_lora_rank=8, lora_dtype=torch.float16
    )
    punica_wrapper = get_punica_wrapper(8192, 256, device, lora_config=lora_config)
    assert check_punica_wrapper(punica_wrapper)

    def create_random_embedding_layer():
        embedding = VocabParallelEmbedding(vocab_size, 256)
        embedding.weight.data = torch.rand_like(embedding.weight.data)
        embedding.weight.data[vocab_size:, :] = 0
        lora_embedding = VocabParallelEmbeddingWithLoRA(embedding)
        lora_embedding.create_lora_weights(max_loras, lora_config)

        return embedding, lora_embedding

    for i in range(NUM_RANDOM_SEEDS):
        set_random_seed(i)

        id_to_index = get_random_id_to_index(num_loras, max_loras)
        embedding, lora_embedding = create_random_embedding_layer()
        lora_embedding.set_mapping(punica_wrapper)
        lora_dict, _ = populate_loras(
            id_to_index,
            layer=lora_embedding,
            layer_weights=embedding.weight.T,
        )

        inputs, index_mapping, prompt_mapping = create_random_inputs(
            active_lora_ids=list(lora_dict.keys()),
            num_inputs=num_loras * 3,
            input_size=(200,),
            input_range=(1, vocab_size),
            device=device,
        )
        lora_mapping = LoRAMapping(index_mapping, prompt_mapping, is_prefill=stage)
        punica_wrapper.update_metadata(
            lora_mapping,
            id_to_index,
            max_loras,
            vocab_size,
        )

        lora_result = lora_embedding(torch.cat(inputs))

        expected_results: list[torch.Tensor] = []
        for input_, lora_id in zip(inputs, prompt_mapping):
            lora = lora_dict[lora_id]
            result = embedding(input_)
            after_a = F.embedding(
                input_,
                lora.lora_a.T,
            )
            result += after_a @ lora.lora_b.T
            expected_results.append(result)
        expected_result = torch.cat(expected_results)

        rtol, atol = TOLERANCES[lora_result.dtype]
        torch.testing.assert_close(lora_result, expected_result, rtol=rtol, atol=atol)

        # Check that resetting the lora weights succeeds

        for slot_idx in range(max_loras):
            lora_embedding.reset_lora(slot_idx)

        inputs, index_mapping, prompt_mapping = create_random_inputs(
            active_lora_ids=[0],
            num_inputs=num_loras * 3,
            input_size=(200,),
            input_range=(1, vocab_size),
            device=device,
        )
        lora_mapping = LoRAMapping(index_mapping, prompt_mapping, is_prefill=stage)
        punica_wrapper.update_metadata(
            lora_mapping,
            id_to_index,
            max_loras,
            vocab_size,
        )

        lora_result = lora_embedding(torch.cat(inputs))
        expected_result = embedding(torch.cat(inputs))

        rtol, atol = TOLERANCES[lora_result.dtype]
        torch.testing.assert_close(lora_result, expected_result, rtol=rtol, atol=atol)