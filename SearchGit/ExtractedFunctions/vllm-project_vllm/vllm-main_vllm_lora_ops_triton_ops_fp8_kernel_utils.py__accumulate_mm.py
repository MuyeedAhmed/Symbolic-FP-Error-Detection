def _accumulate_mm(
    tiled_a,
    tiled_b,
    accumulator,
    a_scale_ptr,
    b_scale_ptr,
    a_scale_k_stride,
    b_scale_k_stride,
    iter_k,
    group_k: tl.constexpr,
    group_n: tl.constexpr,
    use_fp8_w8a8: tl.constexpr,
):
    """
    Core matrix multiplication and accumulation logic with quantization support.

    Args:
        tiled_a (tl.tensor): Loaded tile from A matrix
        tiled_b (tl.tensor): Loaded tile from B matrix
        accumulator (tl.tensor): Current accumulator value
        a_scale_ptr (tl.tensor): Scale pointer for A matrix
        b_scale_ptr (tl.tensor): Scale pointer for B matrix
        a_scale_k_stride (int): K dimension stride for A's block-wise scales
        b_scale_k_stride (int): K dimension stride for B's block-wise scales
        iter_k (int): Current iteration's global K offset
        group_k: Block size for K dimension in block-wise quantization
        group_n: Block size for N dimension in block-wise quantization
        use_fp8_w8a8: Whether using FP8 W8A8 quantization
    """

    if use_fp8_w8a8:
        if group_k > 0 and group_n > 0:
            # Block-wise quantization: scales are loaded per block
            offs_ks = iter_k // group_k
            # a_scale_ptr is (BLOCK_M,) tensor of base pointers per row
            # Load scale for current K-group, result shape: (BLOCK_M,)
            a_scale = tl.load(a_scale_ptr + offs_ks * a_scale_k_stride)
            # b_scale_ptr is (BLOCK_N,) tensor with N-offset pre-baked
            # Load scale for current K-group, result shape: (BLOCK_N,)
            b_scale = tl.load(b_scale_ptr + offs_ks * b_scale_k_stride)
            accumulator += (
                tl.dot(tiled_a, tiled_b) * a_scale[:, None] * b_scale[None, :]
            )
        else:
            # Tensor-wise or per-channel: accumulate and scale at end
            accumulator = tl.dot(tiled_a, tiled_b, acc=accumulator)
    else:
        accumulator += tl.dot(tiled_a, tiled_b)
    return accumulator