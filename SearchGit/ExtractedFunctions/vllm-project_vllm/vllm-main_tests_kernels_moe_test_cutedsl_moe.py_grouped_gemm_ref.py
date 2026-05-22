def grouped_gemm_ref(
    hidden_states_expanded: torch.Tensor,
    hidden_states_3d: torch.Tensor,
    weights: torch.Tensor,
    topk_idx: torch.Tensor,
    masked_m: torch.Tensor,
    B: int,
    topk: int,
    num_experts: int,
    *,
    block_size: int = 16,
) -> torch.Tensor:
    """
    Computes the reference grouped GEMM (fp4 quantized per-expert loop),
    computes flashinfer grouped GEMM (for scale consistency),
    and returns ONLY the repacked reference output: out_ref.

    Returns:
        out_ref: Tensor [num_experts, max_m, n_out]
    """
    device_hs = hidden_states_expanded.device
    device_w = weights.device
    out_dtype = weights.dtype
    n_out = weights.shape[1]

    # Flattened reference output (B*topk, n_out)
    out = torch.zeros((B * topk, n_out), dtype=out_dtype, device=device_w)

    # Per-expert reference compute loop
    for i in range(num_experts):
        mask = topk_idx.view(-1) == i
        if mask.any():
            lhs = hidden_states_expanded[mask]
            rhs = weights[i]

            a_amax = lhs.abs().max().to(torch.float32).to(device_hs)
            b_amax = rhs.abs().max().to(torch.float32).to(device_w)

            a_gs = FLOAT8_E4M3_MAX * FLOAT4_E2M1_MAX / a_amax
            b_gs = FLOAT8_E4M3_MAX * FLOAT4_E2M1_MAX / b_amax

            lhsq, lhsq_sf = fp4_quantize(lhs, a_gs)
            rhsq, rhsq_sf = fp4_quantize(rhs, b_gs)

            lhs_in_dtype = dequantize_nvfp4_to_dtype(
                lhsq,
                lhsq_sf,
                a_gs,
                dtype=lhs.dtype,
                device=device_hs,
                block_size=block_size,
            )
            rhs_in_dtype = dequantize_nvfp4_to_dtype(
                rhsq,
                rhsq_sf,
                b_gs,
                dtype=rhs.dtype,
                device=device_w,
                block_size=block_size,
            )

            out[mask] = lhs_in_dtype @ rhs_in_dtype.t()

    # Determine per-expert max_m
    max_m_val = int(masked_m.max().item())

    # Repack into [num_experts, max_m, n_out]
    out_ref = torch.zeros(
        (num_experts, max_m_val, n_out),
        dtype=out.dtype,
        device=out.device,
    )
    expert_slot = [0] * num_experts

    for i, expert_id in enumerate(topk_idx.view(-1).tolist()):
        slot = expert_slot[expert_id]
        if slot < max_m_val:
            out_ref[expert_id, slot, :] = out[i]
            expert_slot[expert_id] += 1
        else:
            raise IndexError(
                f"Expert {expert_id} exceeded max slots ({max_m_val}). "
                "Increase max_m or check masked_m."
            )

    return out_ref