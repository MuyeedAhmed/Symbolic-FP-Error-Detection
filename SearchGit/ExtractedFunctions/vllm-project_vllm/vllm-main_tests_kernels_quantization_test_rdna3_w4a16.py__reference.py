def _reference(
    x_mk: torch.Tensor,
    q_int4_kn: torch.Tensor,
    scales_gn: torch.Tensor,
    zeros_gn: torch.Tensor | None,
    group_size: int,
    bias: torch.Tensor | None,
) -> torch.Tensor:
    """fp32 reference for the RDNA3 W4A16 op.

    x_mk:       [M, K] fp16/bf16 activations.
    q_int4_kn:  [K, N] int32 raw stored nibbles in [0, 15].
    scales_gn:  [K//G, N] per-group scales (act dtype).
    zeros_gn:   [K//G, N] int32 raw stored zero points in [0, 15], or None
                for the symmetric path (kernel synthesizes stored zero = 7).
    group_size: G.

    The kernel applies the GPTQv1 "+1" zero-point quirk, so the effective
    zero is ``stored_zero + 1`` (symmetric path: 7 + 1 == bias == 8).
    """
    K, N = q_int4_kn.shape
    s_full = scales_gn.repeat_interleave(group_size, dim=0).to(torch.float32)  # [K,N]
    if zeros_gn is None:
        z_full = torch.full(
            (K, N), float(WEIGHT_TYPE.bias), device=x_mk.device, dtype=torch.float32
        )
    else:
        z_full = (zeros_gn + 1).repeat_interleave(group_size, dim=0).to(torch.float32)
    w_fp = (q_int4_kn.to(torch.float32) - z_full) * s_full  # [K, N]
    out = x_mk.to(torch.float32) @ w_fp  # [M, N]
    if bias is not None:
        out = out + bias.to(torch.float32)
    return out.to(x_mk.dtype)