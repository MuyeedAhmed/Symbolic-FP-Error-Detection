def _group_gemm_same_nk_kernel(
    a_ptr,
    b_ptr,
    c_ptr,
    cumsum_M,
    max_M,
    G: tl.constexpr,
    N: tl.constexpr,
    K: tl.constexpr,
    TRANSPOSE_B: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
    BLOCK_K: tl.constexpr,
    GROUP: tl.constexpr,
):
    pid_m, pid_n = _get_pid_mn(tl.program_id(0), max_M, N, BLOCK_M, BLOCK_N, GROUP)
    gid = tl.program_id(1).to(tl.int64)

    gtid_start = tl.load(cumsum_M + gid - 1, mask=gid > 0, other=0).to(tl.int64)
    gtid_end = tl.load(cumsum_M + gid).to(tl.int64)
    m_size = gtid_end - gtid_start

    if pid_m * BLOCK_M >= m_size:
        return

    offs_m = pid_m * BLOCK_M + tl.arange(0, BLOCK_M)
    offs_n = pid_n * BLOCK_N + tl.arange(0, BLOCK_N)
    offs_k = tl.arange(0, BLOCK_K)

    # a is (total_M, K) row-major, offset by expert start
    a_base = a_ptr + gtid_start * K
    # b is (G, N, K) if TRANSPOSE_B else (G, K, N)
    b_base = b_ptr + gid * K * N
    # c is (total_M, N) row-major
    c_base = c_ptr + gtid_start * N

    if TRANSPOSE_B:
        # b layout: (G, N, K), we compute a @ b.T = a(M,K) @ b(N,K).T -> (M,N)
        a_ptrs = a_base + offs_m[:, None] * K + offs_k[None, :]
        b_ptrs = b_base + offs_n[:, None] * K + offs_k[None, :]
    else:
        # b layout: (G, K, N), we compute a @ b = a(M,K) @ b(K,N) -> (M,N)
        a_ptrs = a_base + offs_m[:, None] * K + offs_k[None, :]
        b_ptrs = b_base + offs_k[:, None] * N + offs_n[None, :]

    acc = tl.zeros((BLOCK_M, BLOCK_N), dtype=tl.float32)

    for k_start in range(0, K, BLOCK_K):
        k_offs = k_start + offs_k
        k_mask = k_offs < K

        a_block = tl.load(a_ptrs, mask=(offs_m[:, None] < m_size) & k_mask[None, :], other=0.0)

        if TRANSPOSE_B:
            b_block = tl.load(b_ptrs, mask=(offs_n[:, None] < N) & k_mask[None, :], other=0.0)
            acc += tl.dot(a_block, tl.trans(b_block))
        else:
            b_block = tl.load(b_ptrs, mask=k_mask[:, None] & (offs_n[None, :] < N), other=0.0)
            acc += tl.dot(a_block, b_block)

        if TRANSPOSE_B:
            a_ptrs += BLOCK_K
            b_ptrs += BLOCK_K
        else:
            a_ptrs += BLOCK_K
            b_ptrs += BLOCK_K * N

    c_ptrs = c_base + offs_m[:, None] * N + offs_n[None, :]
    c_mask = (offs_m[:, None] < m_size) & (offs_n[None, :] < N)
    tl.store(c_ptrs, acc.to(c_ptr.dtype.element_ty), mask=c_mask)