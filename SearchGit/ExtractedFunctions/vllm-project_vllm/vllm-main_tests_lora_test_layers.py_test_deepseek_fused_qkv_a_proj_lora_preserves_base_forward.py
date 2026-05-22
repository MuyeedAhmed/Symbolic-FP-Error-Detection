def test_deepseek_fused_qkv_a_proj_lora_preserves_base_forward(
    default_vllm_config, dist_init, device, stage, fully_sharded
):
    if current_platform.is_cuda_alike() or current_platform.is_xpu():
        torch.accelerator.set_device_index(device)

    torch.set_default_device(device)
    dtype = (
        torch.float16
        if (current_platform.is_cuda_alike() or current_platform.is_xpu())
        else torch.float32
    )
    max_loras = 8
    lora_config = LoRAConfig(
        max_loras=max_loras,
        max_lora_rank=8,
        lora_dtype=dtype,
        fully_sharded_loras=fully_sharded,
    )
    punica_wrapper = get_punica_wrapper(8192, 256, device, lora_config=lora_config)
    assert check_punica_wrapper(punica_wrapper)

    class OffsetDeepSeekFusedQkvAProjLinear(DeepSeekV2FusedQkvAProjLinear):
        def forward(self, input_):
            output, output_bias = super().forward(input_)
            return output + 1, output_bias

    layer = OffsetDeepSeekFusedQkvAProjLinear(
        32, [16, 16], prefix="model.layers.0.self_attn.fused_qkv_a_proj"
    )
    layer.weight.data = torch.rand_like(layer.weight.data, dtype=dtype)

    lora_layer = MergedColumnParallelLinearWithLoRA(layer)
    lora_layer.create_lora_weights(max_loras, lora_config)
    lora_layer.set_mapping(punica_wrapper)

    id_to_index = get_random_id_to_index(1, max_loras, log=False)
    active_slot = next(i for i, lora_id in enumerate(id_to_index) if lora_id == 1)
    lora_a = [
        torch.rand(8, 32, dtype=dtype, device=device),
        torch.rand(8, 32, dtype=dtype, device=device),
    ]
    lora_b = [
        torch.rand(16, 8, dtype=dtype, device=device),
        torch.rand(16, 8, dtype=dtype, device=device),
    ]
    lora_layer.set_lora(active_slot, lora_a=lora_a, lora_b=lora_b)

    inputs, index_mapping, prompt_mapping = create_random_inputs(
        active_lora_ids=[1],
        num_inputs=4,
        input_size=(1, 32),
        input_range=(0, 1),
        input_type=dtype,
        device=device,
    )
    lora_mapping = LoRAMapping(index_mapping, prompt_mapping, is_prefill=stage)
    punica_wrapper.update_metadata(lora_mapping, id_to_index, max_loras, 512)

    lora_result = lora_layer(torch.cat(inputs))[0]

    expected_results = []
    for input_ in inputs:
        result = layer(input_)[0]
        result[:, :16] += input_ @ lora_a[0].T @ lora_b[0].T
        result[:, 16:] += input_ @ lora_a[1].T @ lora_b[1].T
        expected_results.append(result)

    rtol, atol = TOLERANCES[lora_result.dtype]
    torch.testing.assert_close(
        lora_result, torch.cat(expected_results), rtol=rtol, atol=atol
    )

    merged_layer = OffsetDeepSeekFusedQkvAProjLinear(
        32, [16, 16], prefix="model.layers.0.self_attn.fused_qkv_a_proj"
    )
    merged_layer.weight.data = layer.weight.data.clone()
    merged_layer.weight.data[:16].add_(lora_b[0] @ lora_a[0])
    merged_layer.weight.data[16:].add_(lora_b[1] @ lora_a[1])
    merged_result = merged_layer(torch.cat(inputs))[0]

    torch.testing.assert_close(lora_result, merged_result, rtol=rtol, atol=atol)