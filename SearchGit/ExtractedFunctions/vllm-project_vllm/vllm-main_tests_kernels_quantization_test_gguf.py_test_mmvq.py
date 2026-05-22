def test_mmvq(hidden_size: int, dtype: torch.dtype, quant_type: GGMLQuantizationType):
    set_random_seed(0)

    tensors = get_gguf_sample_tensors(hidden_size, quant_type)
    x = torch.rand((1, hidden_size), dtype=dtype, device="cuda")
    for tensor in tensors:
        weight = torch.tensor(dequantize(tensor.data, quant_type), device="cuda").to(
            dtype
        )
        ref_output = x @ weight.T

        qweight = torch.tensor(tensor.data, device="cuda")
        output = ops.ggml_mul_mat_vec_a8(qweight, x, quant_type, qweight.shape[0]).to(
            dtype
        )

        torch.testing.assert_close(output, ref_output, atol=1, rtol=1e-1)