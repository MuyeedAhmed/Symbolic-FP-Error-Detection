def test_fused_quantization_quest(rot_size: int):
    dtype, device = DTYPE, DEVICE
    h = get_hadamard_matrix(rot_size, dtype, device)
    x = torch.randn(2, 4096, 4096, dtype=dtype, device=device) * 25.0

    xh_dq_ref, _, _ = _forward_quantize_ref(x, h, rot_size, quest=True)
    xh_e2m1, xh_e8m0 = fusedQuantizeMx(x, h, method="quest")
    xh_e8m0 = xh_e8m0.reshape(2, 4096, 4096 // 32)
    xh_dq, *_ = _dq_fp4(xh_e2m1, xh_e8m0, alpha=1.0)

    torch.testing.assert_close(xh_dq, xh_dq_ref, rtol=0.34, atol=100)
    assert (xh_dq != xh_dq_ref).float().mean() <= 1e-4

    m, n, k = 504, 504, 2048
    a = torch.randn(m, k, dtype=dtype, device=device) * 25.0
    b = torch.randn(n, k, dtype=dtype, device=device) * 25.0

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