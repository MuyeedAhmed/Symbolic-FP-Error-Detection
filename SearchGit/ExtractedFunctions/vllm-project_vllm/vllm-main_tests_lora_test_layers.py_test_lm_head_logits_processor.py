def test_lm_head_logits_processor(
    default_vllm_config, dist_init, num_loras, device, vocab_size, stage
) -> None:
    if current_platform.is_cuda_alike() or current_platform.is_xpu():
        torch.accelerator.set_device_index(device)

    torch.set_default_device(device)
    max_loras = 8
    lora_config = LoRAConfig(
        max_loras=max_loras, max_lora_rank=8, lora_dtype=torch.float16
    )
    punica_wrapper = get_punica_wrapper(8192, 256, device, lora_config=lora_config)
    assert check_punica_wrapper(punica_wrapper)

    def _pretest():
        linear = ParallelLMHead(
            num_embeddings=vocab_size,
            embedding_dim=1024,
            params_dtype=torch.float16,
        )
        linear.weight.data = torch.rand_like(linear.weight.data)
        linear.weight.data[:, vocab_size:] = 0
        logits_processor = LogitsProcessor(vocab_size)
        lora_logits_processor = LogitsProcessorWithLoRA(
            logits_processor, 1024, linear.weight.dtype, linear.weight.device, None
        )
        lora_logits_processor.create_lora_weights(max_loras, lora_config)

        return linear, logits_processor, lora_logits_processor

    for i in range(NUM_RANDOM_SEEDS):
        set_random_seed(i)

        id_to_index = get_random_id_to_index(num_loras, max_loras)
        linear, logits_processor, lora_logits_processor = _pretest()
        lora_logits_processor.set_mapping(punica_wrapper)

        lora_dict, _ = populate_loras(
            id_to_index,
            layer=lora_logits_processor,
            layer_weights=linear.weight,
        )

        inputs, index_mapping, prompt_mapping = create_random_inputs(
            active_lora_ids=list(lora_dict.keys()),
            num_inputs=8 * num_loras,  # * 3,
            input_size=(1, 1024),
            input_range=(0, 1),
            input_type=torch.float16,
            device=device,
        )
        lora_mapping = LoRAMapping(index_mapping, prompt_mapping, is_prefill=stage)
        punica_wrapper.update_metadata(
            lora_mapping,
            id_to_index,
            max_loras,
            vocab_size,
        )
        input_ = torch.rand(20, 1024)

        lora_result = lora_logits_processor._get_logits(
            hidden_states=torch.cat(inputs), lm_head=linear, embedding_bias=None
        )

        original_lm_head = deepcopy(linear)

        expected_results: list[torch.Tensor] = []
        for input_, lora_id in zip(inputs, prompt_mapping):
            lora = lora_dict[lora_id]
            result = logits_processor._get_logits(
                hidden_states=input_, lm_head=linear, embedding_bias=None
            )

            result += input_ @ lora.lora_a.T @ lora.lora_b.T * lora.scaling
            expected_results.append(result)
        expected_result = torch.cat(expected_results)

        # Check that resetting the lora weights succeeds

        for slot_idx in range(max_loras):
            lora_logits_processor.reset_lora(slot_idx)

        inputs, index_mapping, prompt_mapping = create_random_inputs(
            active_lora_ids=[0],
            num_inputs=8 * num_loras * 3,
            input_size=(1, 1024),
            input_range=(0, 1),
            input_type=torch.float16,
            device=device,
        )
        lora_mapping = LoRAMapping(index_mapping, prompt_mapping, is_prefill=stage)
        punica_wrapper.update_metadata(
            lora_mapping,
            id_to_index,
            max_loras,
            vocab_size,
        )

        lora_result = lora_logits_processor._get_logits(
            hidden_states=torch.cat(inputs),
            lm_head=original_lm_head,
            embedding_bias=None,
        )[:, :vocab_size]
        expected_result = logits_processor._get_logits(
            hidden_states=torch.cat(inputs),
            lm_head=original_lm_head,
            embedding_bias=None,
        )

        rtol, atol = TOLERANCES[lora_result.dtype]
        torch.testing.assert_close(lora_result, expected_result, rtol=rtol, atol=atol)