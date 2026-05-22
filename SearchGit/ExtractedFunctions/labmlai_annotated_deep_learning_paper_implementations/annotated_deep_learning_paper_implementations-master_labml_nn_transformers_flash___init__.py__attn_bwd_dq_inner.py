def _attn_bwd_dq_inner(b_dq, b_q, p_kT, p_vT,
                       b_do, b_lse, b_pdp,
                       BLOCK_Q: tl.constexpr, BLOCK_K: tl.constexpr,
                       i, j, steps,
                       MASK: tl.constexpr,
                       q_seq_len: tl.constexpr,
                       kv_seq_len: tl.constexpr):
    """
    #### Inner loop to calculate $dQ_i$
    """

    # Offsets
    offs_i = i + tl.arange(0, BLOCK_Q)
    offs_j = j + tl.arange(0, BLOCK_K)

    # Move the pointers
    p_kT = tl.advance(p_kT, (0, j))
    p_vT = tl.advance(p_vT, (0, j))

    tl.static_assert(BLOCK_Q % BLOCK_K == 0, 'BLOCK_Q must be divisible by BLOCK_K')

    # Iterate over $K$
    for _ in range(steps):
        # Load $K_j^T$
        b_kT = tl.load(p_kT, boundary_check=(1,), padding_option="zero")
        # Load $V_j^T$
        b_vT = tl.load(p_vT, boundary_check=(1,), padding_option="zero")

        # $(\log_2 e) S_{ij} = \sigma (\log_2 e) Q_i K_j^T$
        b_s = tl.dot(b_q, b_kT, out_dtype=HI_PRES_TL)

        # \begin{align}
        # P_{ij} &= \frac{e^{S_{ij}}}{L_i}
        # \\
        # &= \frac{2^{(log_2 e) S_{ij}}}{2^{\log_2 L_i}}
        # \\
        # &= 2^{(log_2 e) S_{ij} - \log_2 L_i}
        # \end{align}
        b_p = tl.math.exp2(b_s - b_lse[:, None])

        # Autoregressive masking
        if MASK:
            causal_mask = (offs_i[:, None] >= offs_j[None, :])
            b_p = tl.where(causal_mask, b_p, 0.0)

        # Mask out if the block is beyond the end of $Q_i$
        j_mask = offs_j < kv_seq_len
        b_p = tl.where(j_mask[None, :], b_p, 0.0)

        # $$dq_i = \sum_j dS_{ij} k_j = \sum_j P_{ij} \big( dP_{ij} - D_i \big) k_j$$

        # $dP_{ij} = dO_i V_j^T$
        b_dp = tl.dot(b_do, b_vT, out_dtype=HI_PRES_TL).to(HI_PRES_TL)
        # $dS_{ij} = P_{ij} \big( dP_{ij} - D_i \big)$
        b_ds = b_p * (b_dp - b_pdp[:, None])
        # $(\log_2 e) dQ_i = \sum_j dS_{ij} \sigma (\log_2 e) K_j$
        b_dq += tl.dot(b_ds.to(b_kT.dtype), tl.trans(b_kT), out_dtype=HI_PRES_TL)

        # Increment pointers.
        offs_j += BLOCK_K
        p_kT = tl.advance(p_kT, (0, BLOCK_K))
        p_vT = tl.advance(p_vT, (0, BLOCK_K))

    # Return accumulated $dQ$
    return b_dq