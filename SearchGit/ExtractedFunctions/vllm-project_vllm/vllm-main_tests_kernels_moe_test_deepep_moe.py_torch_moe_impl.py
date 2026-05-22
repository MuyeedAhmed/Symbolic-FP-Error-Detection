def torch_moe_impl(
    test_tensors: TestTensors,
    w1: torch.Tensor,
    w2: torch.Tensor,
    w1_scale: torch.Tensor | None,
    w2_scale: torch.Tensor | None,
    using_fp8_dispatch: bool,
    per_act_token_quant: bool,
):
    a, topk_ids, topk_weights = (
        test_tensors.rank_tokens,
        test_tensors.topk,
        test_tensors.topk_weights,
    )
    if using_fp8_dispatch:
        # The DeepEP implementation is requested to dispatch using FP8.
        # For numerical stability for testing, emulate the fp8 dispatch by
        # blockwise quant and de-quant.
        assert not per_act_token_quant
        a = test_tensors.rank_tokens
        aq, aq_scale = per_token_group_quant_fp8(a, 128, use_ue8m0=False)
        a = (
            (aq.view(-1, 128).to(torch.float32) * aq_scale.view(-1, 1))
            .view(a.shape)
            .to(a.dtype)
        )

    is_quantized = w1.dtype == torch.float8_e4m3fn
    a_dtype = a.dtype
    if is_quantized:
        w1 = w1.to(dtype=torch.float32) * w1_scale
        w2 = w2.to(dtype=torch.float32) * w2_scale
        a = a.to(dtype=torch.float32)

    m, _ = a.shape
    topk = topk_ids.size(1)
    out = torch.zeros_like(a)

    for i in range(m):
        a_i = a[i]
        o_i = out[i]
        for j in range(topk):
            e = topk_ids[i][j]
            e_w = topk_weights[i][j]
            w1_e = w1[e]
            w2_e = w2[e]
            o_i += (
                SiluAndMul()(a_i @ w1_e.transpose(0, 1)) @ w2_e.transpose(0, 1)
            ) * e_w

    if is_quantized:
        out = out.to(dtype=a_dtype)

    return out