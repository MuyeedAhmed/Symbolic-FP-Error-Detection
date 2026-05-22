def _attn_bwd_dkdv_inner(b_dk, b_dv,
                         p_qT, b_k, b_v, p_do,
                         p_lse, p_pdp,
                         BLOCK_Q: tl.constexpr, BLOCK_K: tl.constexpr,
                         d_head: tl.constexpr,
                         j, i, steps,
                         MASK: tl.constexpr,
                         q_seq_len: tl.constexpr,
                         kv_seq_len: tl.constexpr):
    """
    #### Inner loop to calculate $dK_j$, $dV_j$
    """

    # To apply the mask
    tl.static_assert(BLOCK_K % BLOCK_Q == 0)

    # Offsets and mask
    offs_i = i + tl.arange(0, BLOCK_Q)
    offs_j = j + tl.arange(0, BLOCK_K)

    # Move the pointers
    p_qT = tl.advance(p_qT, (0, i))
    p_do = tl.advance(p_do, (i, 0))
    p_lse = tl.advance(p_lse, (i,))
    p_pdp = tl.advance(p_pdp, (i,))

    # Iterate over $Q$
    for _ in range(steps):
        # Load $Q_i^T$
        b_qT = tl.load(p_qT, boundary_check=(1,), padding_option="zero")

        # $log_2 L_i$
        b_l = tl.load(p_lse, boundary_check=(0,), padding_option="zero")

        # $(\log_2 e) S_{ij}^T = \sigma (\log_2 e) K_j Q_i^T$
        b_sT = tl.dot(b_k, b_qT, out_dtype=HI_PRES_TL)

        # \begin{align}
        # P_{ij} &= \frac{e^{S_{ij}}}{L_i}
        # \\
        # &= \frac{2^{(log_2 e) S_{ij}}}{2^{\log_2 L_i}}
        # \\
        # &= 2^{(log_2 e) S_{ij} - \log_2 L_i}
        # \end{align}
        b_pT = tl.math.exp2(b_sT - b_l[None, :])

        # Autoregressive masking
        if MASK:
            mask = (offs_i[None, :] >= offs_j[:, None])
            b_pT = tl.where(mask, b_pT, 0.0)

        # Mask out if the block is beyond the end of $Q_i$
        #
        # Note: No need to mask out based on $j$
        # because the effects on positions outside boundary will not get stored in $dK$ or $dV$
        # Masking by $i$ may also not be necessary size the tensors have 0 on loading
        i_mask = offs_i < q_seq_len
        b_pT = tl.where(i_mask[None, :], b_pT, 0.0)

        # $dV_j = \sum_i P_{ij} dO_i$
        b_do = tl.load(p_do, boundary_check=(0,), padding_option="zero")
        b_dv += tl.dot(b_pT.to(b_do.dtype), b_do, out_dtype=HI_PRES_TL)

        # $D_i$
        b_pdp = tl.load(p_pdp, boundary_check=(0,), padding_option="zero")
        # $dP_{ij} = V_j dO_i^T$
        b_dpT = tl.dot(b_v, tl.trans(b_do), out_dtype=HI_PRES_TL).to(HI_PRES_TL)
        # $dS_{ij} = P_{ij} \big( dP_{ij} - D_i \big)$
        b_dsT = b_pT * (b_dpT - b_pdp[None, :])
        # $\frac{1}{\sigma} dK_j = \sum_i dS_{ij} Q_i$
        b_dk += tl.dot(b_dsT.to(b_qT.dtype), tl.trans(b_qT), out_dtype=HI_PRES_TL)

        # Increment pointers.
        offs_i += BLOCK_Q
        p_lse = tl.advance(p_lse, (BLOCK_Q,))
        p_pdp = tl.advance(p_pdp, (BLOCK_Q,))
        p_qT = tl.advance(p_qT, (0, BLOCK_Q))
        p_do = tl.advance(p_do, (BLOCK_Q, 0))

    # Return accumulated $dK$ and $dV$
    return b_dk, b_dv