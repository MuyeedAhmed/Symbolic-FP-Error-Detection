def test_cutlass_int8_azp_bias_fold(m: int, n: int, k: int, out_dtype: torch.dtype):
    # Currently, the test is failing because folding azp into
    # 16-bit bias loses too much precision
    scale_a = torch.randn((1, 1), device="cuda", dtype=torch.float32) / 10
    scale_b = torch.randn((1, n), device="cuda", dtype=torch.float32) / 10

    aq_i8 = rand_int8((m, k))
    bq_i8 = rand_int8((n, k)).t()

    aq_i32 = aq_i8.to(dtype=torch.int32)
    bq_i32 = bq_i8.to(dtype=torch.int32)

    aq_f32 = aq_i8.to(dtype=torch.float32)
    bq_f32 = bq_i8.to(dtype=torch.float32)

    b_dq = scale_b * bq_f32

    azp_a = torch.rand((1,), device="cuda", dtype=torch.float32) * 10 + 1.5
    azp_aq_i8 = (azp_a / scale_a).to(dtype=torch.int8)
    azp_a = azp_aq_i8.to(dtype=torch.float32) * scale_a  # correct for rounding

    a_dq = scale_a * (aq_i32 + azp_aq_i8).to(dtype=torch.float32)
    torch.testing.assert_close(a_dq, scale_a * aq_f32 + azp_a)

    baseline_dq = torch.mm(a_dq, b_dq).to(out_dtype)

    J = torch.ones((1, k), device="cuda", dtype=torch.float32)
    azp_bias = (azp_a * scale_b * (J @ bq_f32)).to(out_dtype)
    assert azp_bias.shape == (1, n)
    assert azp_bias[0, :].shape == (n,)

    baseline_q = (
        scale_a.to(device="cpu")
        * scale_b.to(device="cpu")
        * ((aq_i32 + azp_aq_i8).to(device="cpu") @ bq_i32.to(device="cpu"))
    ).to(dtype=out_dtype, device="cuda")

    out = ops.cutlass_scaled_mm(
        aq_i8, bq_i8, scale_a, scale_b, out_dtype=out_dtype, bias=azp_bias[0, :]
    )
    torch.testing.assert_close(out, baseline_dq, rtol=1e-2, atol=1e0)
    torch.testing.assert_close(out, baseline_q, rtol=1e-2, atol=1e0)