def _compute_compressed_swizzled_bitmask(dense):
    """
    Calculates the compressed swizzled bitmask from a dense tensor
    """

    # first we need to convert the dense tensor to a bitmask
    int_bitmask = dense.bool().to(torch.uint8)

    # Each thread is responsible for an 8x8 tile, which contains 4 4x4 tiles:
    # A, B, C and D, as displayed in the following schema:
    # +---+---+
    # | A | B |
    # +---+---+
    # | C | D |
    # +---+---+

    # we first need to split into the 8x8 tiles
    bitmask_8x8_chunks = int_bitmask.unfold(0, 8, 8).unfold(1, 8, 8)

    # then we unfold again to get our individual 4x4 tiles
    bitmask_4x4_chunks = bitmask_8x8_chunks.unfold(2, 4, 4).unfold(3, 4, 4)

    # Each 4x4 bitmask defines two 8-bit integers, which encode the sparsity pattern
    # of that tile. Note that the least significant bit is stored first.
    # [1 1 0 0]
    # [1 1 0 0]  ->  0011 0011 ->   51
    # [0 0 1 1]      1100 1100      204
    # [0 0 1 1]

    # reshape tensor to expand tiles into 8-bit vectors
    bitmask_binary_representation = bitmask_4x4_chunks.reshape(
        *bitmask_4x4_chunks.shape[:2], 4, 2, 8
    )

    # to convert from binary representation, we can do a matmul with powers of two
    powers_of_two = 2 ** torch.arange(8, dtype=torch.float, device="cuda")
    # To run on GPU: cast to float to do matmul and then cast back
    compressed_swizzled_bitmask = (
        bitmask_binary_representation.to(torch.float) @ powers_of_two
    ).to(torch.uint8)

    return compressed_swizzled_bitmask