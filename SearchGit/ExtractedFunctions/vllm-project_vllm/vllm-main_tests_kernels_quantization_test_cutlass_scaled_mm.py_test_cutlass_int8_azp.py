def test_cutlass_int8_azp(
    m: int, n: int, k: int, out_dtype: torch.dtype, use_bias: bool, azp_per_token: bool
):
    m_azp = m if azp_per_token else 1
    scale_a = torch.randn((m_azp, 1), device="cuda", dtype=torch.float32) / 10
    scale_b = torch.randn((1, n), device="cuda", dtype=torch.float32) / 10

    aq_i8 = rand_int8((m, k))
    aq_i32 = aq_i8.to(dtype=torch.int32)
    aq_f32 = aq_i8.to(dtype=torch.float32)

    bq_i8 = rand_int8((n, k)).t()
    bq_i32 = bq_i8.to(dtype=torch.int32)
    bq_f32 = bq_i8.to(dtype=torch.float32)
    b_dq = scale_b * bq_f32

    azp_a = torch.rand((m_azp, 1), device="cuda", dtype=torch.float32) * 10 + 1.5
    azp_aq_i8 = (azp_a / scale_a).to(dtype=torch.int8)
    azp_a = azp_aq_i8.to(dtype=torch.float32) * scale_a  # correct for rounding

    a_dq = scale_a * (aq_i32 - azp_aq_i8).to(dtype=torch.float32)
    torch.testing.assert_close(a_dq, scale_a * aq_f32 - azp_a, rtol=1e-4, atol=1e-3)

    if use_bias:
        bias = torch.rand((1, n), device="cuda", dtype=out_dtype) * 10 + 2.5
    else:
        bias = torch.zeros((1, n), device="cuda", dtype=out_dtype)

    baseline_dq = (torch.mm(a_dq, b_dq) + bias).to(out_dtype)

    # int32 mm not supported on CUDA
    a_noazp_i32_cpu = (aq_i32 - azp_aq_i8).to(device="cpu")
    cq = (a_noazp_i32_cpu @ bq_i32.to(device="cpu")).to(device="cuda")
    baseline_q = (scale_a * scale_b * cq + bias).to(dtype=out_dtype)

    # Hadamard is just the sum of the cols
    azp_adj_i32 = bq_i32.sum(dim=0, keepdim=True, dtype=torch.int32)
    azp_i32 = azp_aq_i8.to(dtype=torch.int32)
    func_bias = bias if use_bias else None

    if azp_per_token:
        out = ops.cutlass_scaled_mm_azp(
            aq_i8, bq_i8, scale_a, scale_b, out_dtype, azp_adj_i32, azp_i32, func_bias
        )
    else:
        azp_with_adj_i32 = azp_i32 * azp_adj_i32
        out = ops.cutlass_scaled_mm_azp(
            aq_i8, bq_i8, scale_a, scale_b, out_dtype, azp_with_adj_i32, None, func_bias
        )

    # bfloat16 precision is 7-bit mantissa -> 2^-8 ~ 0.4%
    # float16 precision is 10-bit mantissa -> 2^-11 ~ 0.05%
    rtol = 1e-2 if out_dtype == torch.bfloat16 else 1e-3
    atol = 1e-3
    torch.testing.assert_close(out, baseline_dq, rtol=rtol, atol=atol)
    torch.testing.assert_close(out, baseline_q, rtol=rtol, atol=atol)

    if azp_per_token:
        opcheck(
            torch.ops._C.cutlass_scaled_mm_azp,
            (out, aq_i8, bq_i8, scale_a, scale_b, azp_adj_i32, azp_i32, func_bias),
        )
    else:
        opcheck(
            torch.ops._C.cutlass_scaled_mm_azp,
            (out, aq_i8, bq_i8, scale_a, scale_b, azp_with_adj_i32, None, func_bias),
        )