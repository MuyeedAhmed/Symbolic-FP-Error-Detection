def test_llama_shapes(model: str, layer_idx: int, batch: int, had_size: int):
    dtype, device = DTYPE, DEVICE
    m = batch
    k, n = LLAMA_MODELS[model][layer_idx]

    h = get_hadamard_matrix(had_size, dtype, device)

    a = torch.rand(m, k, dtype=dtype, device=device) * 25.0
    b = torch.rand(n, k, dtype=dtype, device=device) * 25.0

    a_e2m1, a_e8m0 = fusedQuantizeMx(a, h, method="quest")
    b_e2m1, b_e8m0 = fusedQuantizeMx(b, h, method="quest")

    a_dq, *_ = _dq_fp4(a_e2m1, a_e8m0[:m, :k], alpha=1.0)
    b_dq, *_ = _dq_fp4(b_e2m1, b_e8m0[:n, :k], alpha=1.0)
    out_ref = a_dq @ b_dq.transpose(-2, -1)

    a_scale_block = to_blocked(a_e8m0, backend="triton")
    b_scale_block = to_blocked(b_e8m0, backend="triton")
    alpha = torch.tensor([1.0], device=device)
    out = matmul_mxf4_bf16_tn(a_e2m1, b_e2m1, a_scale_block, b_scale_block, alpha)
    assert out.equal(out_ref.to(dtype=out.dtype))