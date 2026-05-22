def _attn_fwd_inner(b_o, b_l, b_m, b_q,
                    p_kT, p_v,
                    sm_scale_log2e,
                    BLOCK_Q: tl.constexpr,
                    d_head: tl.constexpr,
                    BLOCK_K: tl.constexpr,
                    offs_i, offs_j,
                    j,
                    steps,
                    MASK: tl.constexpr,
                    q_seq_len: tl.constexpr,
                    kv_seq_len: tl.constexpr
                    ):
    """
    #### Inner loop to calculate $O_i$

    This iterates through keys and values starting from `j` for `steps` number of steps.
    In each step it processes `BLOCK_K` entries of keys/values.
    """
    tl.static_assert(BLOCK_Q % BLOCK_K == 0)

    # Move $K_j$ and $V_j$ pointers
    p_kT = tl.advance(p_kT, (0, j))
    p_v = tl.advance(p_v, (j, 0))

    # Iterate over $K$, $V$ and update $\tilde{O}_i$ and $l_i$
    for _ in range(steps):
        # Load $K_j^T$
        b_kT = tl.load(p_kT, boundary_check=(1,), padding_option="zero")
        # Compute $(\log 2) S_ij  = (\log 2) \sigma Q_i K_j^T$
        b_s = tl.dot(b_q, b_kT, out_dtype=HI_PRES_TL)
        b_s = b_s * sm_scale_log2e

        # Apply causal mask
        if MASK:
            causal_mask = offs_i[:, None] >= (j + offs_j[None, :])
            b_s = tl.where(causal_mask, b_s, -float("inf"))

        # Mask out if the block is beyond the end of $K_j$
        j_mask = (j + offs_j) < kv_seq_len
        b_s = tl.where(j_mask[None, :], b_s, -float("inf"))

        # $(\log_2 e) m_{i}^{\text{new}} = \max((\log_2 e) m_i, \max_{j=j1}^{j2} (\log_2 e) S_{ij})$
        b_m_new = tl.maximum(b_m, tl.max(b_s, -1))
        # \begin{align}
        # \tilde{P}_{ij} &= e^{(S_{ij} - m_i^{\text{new}}}
        # \\
        # &= 2^{(\log_2 e) S_{ij} - (\log_2 e) m_i^{\text{new}}}
        # \end{align}
        b_p = tl.math.exp2(b_s - b_m_new[:, None])

        # $\sum_{j=j1}^{j2} \tilde{P}_{ij}$
        b_l_new = tl.sum(b_p, -1)
        # $e^{m_i - m_{i}^{\text{new}}}$
        b_m_m_new = tl.math.exp2(b_m - b_m_new)
        # $l_i \leftarrow e^{m_i - m_{i}^{\text{new}}} l_i + \sum_{j=j1}^{j2} \tilde{P}_{ij}$
        b_l = b_l * b_m_m_new + b_l_new

        # $O_i \leftarrow e^{m_i - m_{i}^{\text{new}}} O_i + \tilde{P}_{ij} V_j$
        b_o = b_o * b_m_m_new[:, None]
        b_p = b_p.to(b_q.dtype)  # TODO
        b_v = tl.load(p_v, boundary_check=(0,), padding_option="zero")
        b_o += tl.dot(b_p, b_v, out_dtype=HI_PRES_TL)

        # $(\log_2 e) m_i \leftarrow (\log_2 e) m_{i}^{\text{new}}$
        b_m = b_m_new

        # Move pointers
        j += BLOCK_K
        p_v = tl.advance(p_v, (BLOCK_K, 0))
        p_kT = tl.advance(p_kT, (0, BLOCK_K))

    tl.static_assert(b_o.dtype == HI_PRES_TL, "attn_fwd_inner requires accumulator to be in HI_PRES_TL precision")

    return b_o, b_l, b_m