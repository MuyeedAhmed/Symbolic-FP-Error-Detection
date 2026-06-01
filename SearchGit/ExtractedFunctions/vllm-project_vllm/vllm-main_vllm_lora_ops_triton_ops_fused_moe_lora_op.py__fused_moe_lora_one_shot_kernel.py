def _fused_moe_lora_one_shot_kernel(
    # ---- pointers ----
    x_ptr,
    A_ptrs,
    B_ptrs,
    out_ptr,
    topk_weights_ptr,
    sorted_token_ids_ptr,
    expert_ids_ptr,
    num_tokens_post_padded_ptr,
    token_lora_mapping_ptr,
    lora_ids_ptr,
    adapter_enabled_ptr,
    # ---- dims ----
    N,
    K,
    num_valid_tokens,
    top_k_num,
    max_loras,
    # ---- strides ----
    stride_xm,
    stride_xk,
    stride_A_lora,
    stride_A_expert,
    stride_A_r,
    stride_A_k,
    stride_B_lora,
    stride_B_expert,
    stride_B_n,
    stride_B_r,
    stride_om,
    stride_on,
    stride_tl_,
    stride_el,
    # ---- scalar ----
    slice_n_offset,
    # ---- constexpr (set per call) ----
    token_mapping_factor: tl.constexpr,
    naive_block_assignment: tl.constexpr,
    MUL_ROUTED_WEIGHT: tl.constexpr,
    BLOCK_M: tl.constexpr,
    BLOCK_R: tl.constexpr,
    actual_rank: tl.constexpr,
    NPID_FACTOR: tl.constexpr,
    BLOCK_N: tl.constexpr,
    BLOCK_K: tl.constexpr,
    EVEN_K: tl.constexpr,
    ADD_INPUTS: tl.constexpr,
):
    pid_full = tl.program_id(axis=0)
    pid_m = pid_full // NPID_FACTOR
    pid_n_outer = pid_full % NPID_FACTOR
    slice_id = tl.program_id(axis=1)
    lora_idx = tl.program_id(axis=2)

    # Resolve lora_id.
    if naive_block_assignment:
        token_idx_for_lora = pid_m // top_k_num
        lora_id = tl.load(token_lora_mapping_ptr + token_idx_for_lora)
    else:
        lora_id = tl.load(lora_ids_ptr + lora_idx)
    if lora_id < 0:
        return
    if lora_id >= max_loras:
        return
    enabled = tl.load(adapter_enabled_ptr + lora_id)
    if enabled == 0:
        return

    if not naive_block_assignment:
        ntpp = tl.load(num_tokens_post_padded_ptr + lora_id)
        if pid_m * BLOCK_M >= ntpp:
            return

    # Resolve expert_id.
    if naive_block_assignment:
        expert_id = tl.load(expert_ids_ptr + pid_m)
    else:
        ind = lora_id * stride_el + pid_m
        expert_id = tl.load(
            expert_ids_ptr + ind, mask=ind < max_loras * stride_el, other=-1
        )
    if expert_id < 0:
        return

    # Compute offs_token (flat token ids).
    offs = tl.arange(0, BLOCK_M).to(tl.int64)
    if naive_block_assignment:
        offs_token = tl.where(offs == 0, pid_m, num_valid_tokens)
    else:
        offs_token_id = pid_m * BLOCK_M + offs
        token_ind = stride_tl_ * lora_id + offs_token_id
        offs_token = tl.load(
            sorted_token_ids_ptr + token_ind,
            mask=token_ind < max_loras * stride_tl_,
            other=num_valid_tokens,
        )
    token_mask = offs_token < num_valid_tokens

    # N range owned by this program. Splitting [0, N) into NPID_FACTOR
    # contiguous outer blocks lets us scale parallelism for small batches.
    n_per_outer = tl.cdiv(N, NPID_FACTOR)
    n_lo = pid_n_outer * n_per_outer
    n_hi = tl.minimum((pid_n_outer + 1) * n_per_outer, N)
    if n_lo >= N:
        return

    # Slice pointers.
    cur_A_ptr = tl.load(A_ptrs + slice_id).to(tl.pointer_type(out_ptr.dtype.element_ty))
    cur_B_ptr = tl.load(B_ptrs + slice_id).to(tl.pointer_type(out_ptr.dtype.element_ty))

    A_base = cur_A_ptr + lora_id * stride_A_lora + expert_id * stride_A_expert
    B_base = cur_B_ptr + lora_id * stride_B_lora + expert_id * stride_B_expert

    # SHRINK: tmp[BLOCK_M, BLOCK_R] = x @ A^T, accumulated in fp32 registers.
    offs_r = tl.arange(0, BLOCK_R)
    rank_mask = offs_r < actual_rank
    # Clamp rank offsets so OOB rows of A / B map to address 0; the mask
    # zeros the loaded values. Required when BLOCK_R > actual_rank
    # (e.g. rank=4 padded to 16) -- without clamping, tl.load would address
    # the next expert's memory.
    safe_offs_r = tl.where(rank_mask, offs_r, 0)
    offs_k = tl.arange(0, BLOCK_K)

    offs_x_row = offs_token // token_mapping_factor
    x_ptrs = x_ptr + offs_x_row[:, None] * stride_xm + offs_k[None, :] * stride_xk
    a_ptrs = A_base + offs_k[:, None] * stride_A_k + safe_offs_r[None, :] * stride_A_r

    tmp = tl.zeros((BLOCK_M, BLOCK_R), dtype=tl.float32)
    if EVEN_K:
        for _ in range(0, K, BLOCK_K):
            x = tl.load(x_ptrs, mask=token_mask[:, None], other=0.0)
            a = tl.load(a_ptrs, mask=rank_mask[None, :], other=0.0)
            tmp += tl.dot(x, a)
            x_ptrs += BLOCK_K * stride_xk
            a_ptrs += BLOCK_K * stride_A_k
    else:
        for kb in range(0, K, BLOCK_K):
            k_remain = K - kb
            k_mask = offs_k < k_remain
            x = tl.load(x_ptrs, mask=token_mask[:, None] & k_mask[None, :], other=0.0)
            a = tl.load(a_ptrs, mask=k_mask[:, None] & rank_mask[None, :], other=0.0)
            tmp += tl.dot(x, a)
            x_ptrs += BLOCK_K * stride_xk
            a_ptrs += BLOCK_K * stride_A_k

    tmp_typed = tmp.to(out_ptr.dtype.element_ty)

    # EXPAND: out[tokens, n] += tmp @ B^T, looped over BLOCK_N tiles within
    # this program's [n_lo, n_hi). The (offs_n < n_hi) mask is required
    # whenever BLOCK_N > n_per_outer to keep adjacent outer blocks from
    # writing into each other's columns.
    if MUL_ROUTED_WEIGHT:
        moe_w = tl.load(topk_weights_ptr + offs_token, mask=token_mask, other=0.0).to(
            tl.float32
        )

    out_slice_base = out_ptr + slice_id * slice_n_offset

    for n_start in range(n_lo, n_hi, BLOCK_N):
        offs_n = n_start + tl.arange(0, BLOCK_N)
        n_mask = (offs_n < N) & (offs_n < n_hi)

        b_ptrs = (
            B_base + safe_offs_r[:, None] * stride_B_r + offs_n[None, :] * stride_B_n
        )
        b = tl.load(b_ptrs, mask=rank_mask[:, None] & n_mask[None, :], other=0.0)

        acc = tl.dot(tmp_typed, b)  # (BLOCK_M, BLOCK_N) fp32
        if MUL_ROUTED_WEIGHT:
            acc = acc * moe_w[:, None]

        out_ptrs = (
            out_slice_base
            + offs_token[:, None] * stride_om
            + offs_n[None, :] * stride_on
        )
        out_mask = token_mask[:, None] & n_mask[None, :]
        if ADD_INPUTS:
            prev = tl.load(out_ptrs, mask=out_mask, other=0.0)
            tl.store(out_ptrs, prev + acc.to(out_ptr.dtype.element_ty), mask=out_mask)
        else:
            tl.store(out_ptrs, acc.to(out_ptr.dtype.element_ty), mask=out_mask)