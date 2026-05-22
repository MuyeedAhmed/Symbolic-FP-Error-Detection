def _group_gemm_same_mn_kernel(
    a_ptr,
    b_ptr,
    c_ptr,
    cumsum_K,
    G: tl.constexpr,
    M: tl.constexpr,
    N: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_N: tl.constexpr,
    BLOCK_K: tl.constexpr,
    GROUP: tl.constexpr,
):
    pid_m, pid_n = _get_pid_mn(tl.program_id(0), M, N, BLOCK_M, BLOCK_N, GROUP)
    gid = tl.program_id(1).to(tl.int64)

    gtid_start = tl.load(cumsum_K + gid - 1, mask=gid > 0, other=0).to(tl.int64)
    gtid_end = tl.load(cumsum_K + gid).to(tl.int64)
    k_size = gtid_end - gtid_start

    offs_m = pid_m * BLOCK_M + tl.arange(0, BLOCK_M)
    offs_n = pid_n * BLOCK_N + tl.arange(0, BLOCK_N)

    # c is (G, M, N)
    c_base = c_ptr + gid * M * N

    if k_size == 0:
        c_ptrs = c_base + offs_m[:, None] * N + offs_n[None, :]
        c_mask = (offs_m[:, None] < M) & (offs_n[None, :] < N)
        tl.store(c_ptrs, tl.zeros((BLOCK_M, BLOCK_N), dtype=c_ptr.dtype.element_ty), mask=c_mask)
        return

    acc = tl.zeros((BLOCK_M, BLOCK_N), dtype=tl.float32)
    offs_k = tl.arange(0, BLOCK_K)

    # a is (total_K, M), compute a.T @ b -> (M, N)
    # b is (total_K, N)
    a_base = a_ptr + gtid_start * M
    b_base = b_ptr + gtid_start * N

    for k_start in range(0, k_size, BLOCK_K):
        k_offs = k_start + offs_k
        k_mask = k_offs < k_size

        a_ptrs = a_base + k_offs[:, None] * M + offs_m[None, :]
        a_block_t = tl.trans(tl.load(a_ptrs, mask=k_mask[:, None] & (offs_m[None, :] < M), other=0.0))

        # Load b block: (BLOCK_K, BLOCK_N)
        b_ptrs = b_base + k_offs[:, None] * N + offs_n[None, :]
        b_block = tl.load(b_ptrs, mask=k_mask[:, None] & (offs_n[None, :] < N), other=0.0)

        acc += tl.dot(a_block_t, b_block)

    c_ptrs = c_base + offs_m[:, None] * N + offs_n[None, :]
    c_mask = (offs_m[:, None] < M) & (offs_n[None, :] < N)
    tl.store(c_ptrs, acc.to(c_ptr.dtype.element_ty), mask=c_mask)