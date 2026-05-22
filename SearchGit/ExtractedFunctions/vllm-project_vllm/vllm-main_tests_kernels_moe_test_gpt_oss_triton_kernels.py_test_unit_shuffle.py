def test_unit_shuffle():
    N = ModelConfig.intermediate_size
    K = ModelConfig.hidden_size
    m = torch.randn((K, 2 * N), dtype=torch.bfloat16, device="cuda")

    x = torch.randn(K, dtype=torch.bfloat16, device="cuda")

    m_shuffled = shuffle_weight(m)

    out_ref = x @ m
    out_ref = swiglu(out_ref, limit=1.0)

    out = x @ m_shuffled
    out = triton_kernels.swiglu.swiglu_torch(
        out,
        alpha=1.702,
        precision_config=triton_kernels.swiglu.PrecisionConfig(limit=1.0),
    )

    assert_close(ref=out_ref, tri=out)