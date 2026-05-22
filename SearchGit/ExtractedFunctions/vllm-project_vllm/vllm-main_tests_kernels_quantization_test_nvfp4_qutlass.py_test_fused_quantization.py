def test_fused_quantization(rot_size: int, global_scale_value: float):
    dtype, device = DTYPE, DEVICE
    h = get_hadamard_matrix(rot_size, dtype, device)
    x = torch.randn(2, 4096, 4096, dtype=dtype, device=device) * 25.0
    global_scale = torch.tensor([global_scale_value], device=device)

    xh_dq_ref, _, _ = _forward_quantize_ref(x, h, rot_size)
    xh_e2m1, xh_e4m3 = fusedQuantizeNv(x, h, global_scale)
    xh_e4m3 = xh_e4m3.reshape(2, 4096, 4096 // 16)
    xh_dq, *_ = _dq_fp4(xh_e2m1, xh_e4m3, alpha=global_scale_value)

    torch.testing.assert_close(xh_dq, xh_dq_ref, rtol=0.34, atol=100)
    assert (xh_dq != xh_dq_ref).float().mean() <= 1e-1

    m, n, k = 504, 4096 * 2, 4096
    a = torch.randn(m, k, dtype=dtype, device=device) * 25.0
    b = torch.randn(n, k, dtype=dtype, device=device) * 25.0

    a_e2m1, a_e4m3 = fusedQuantizeNv(a, h, global_scale)
    b_e2m1, b_e4m3 = fusedQuantizeNv(b, h, global_scale)

    a_dq, *_ = _dq_fp4(a_e2m1, a_e4m3[:m, :k], alpha=1.0)
    b_dq, *_ = _dq_fp4(b_e2m1, b_e4m3[:n, :k], alpha=1.0)
    out_ref = a_dq @ b_dq.transpose(-2, -1)

    a_scale_block = to_blocked(a_e4m3, backend="triton").view(-1, k // 16)
    b_scale_block = to_blocked(b_e4m3, backend="triton").view(-1, k // 16)
    alpha = torch.tensor([1.0], device=device)
    out = ops.cutlass_scaled_fp4_mm(
        a_e2m1, b_e2m1, a_scale_block, b_scale_block, alpha, torch.bfloat16
    )
    assert out.equal(out_ref.to(dtype=out.dtype))