def test_llama_shapes(model: str, layer_idx: int, batch: int, rot_size: int):
    dtype, device = DTYPE, DEVICE
    m = batch
    k, n = LLAMA_MODELS[model][layer_idx]

    h = get_hadamard_matrix(rot_size, dtype, device)

    a = torch.randn(m, k, dtype=dtype, device=device) * 25.0
    b = torch.randn(n, k, dtype=dtype, device=device) * 25.0

    global_scale = torch.tensor([1.0], device=device)

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