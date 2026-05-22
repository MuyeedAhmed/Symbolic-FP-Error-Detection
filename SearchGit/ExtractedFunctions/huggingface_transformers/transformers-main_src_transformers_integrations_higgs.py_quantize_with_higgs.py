def quantize_with_higgs(weight, bits: int = 4, p: int = 2, group_size: int = 256, hadamard_size: int = 1024):
    assert len(weight.shape) == 2, "Only 2D weights are supported for now"

    grid = get_higgs_grid(p, 2 ** (p * bits)).to(weight.device)
    grid_norm_2 = torch.linalg.norm(grid, axis=-1) ** 2

    device = weight.device
    dtype = weight.dtype
    weight = weight.to(copy=True, dtype=torch.float32)
    # Pad to Hadamard transform size
    weight = pad_to_block(weight, [1], hadamard_size)

    # Scale and Hadamard transform
    mult = weight.shape[1] // hadamard_size
    weight = weight.reshape(-1, mult, hadamard_size)
    scales = torch.linalg.norm(weight, axis=-1)
    weight = hadamard_transform(weight, 1) / scales[:, :, None]

    # Pad to edenn_d and project
    weight = pad_to_block(weight, [2], p).reshape(weight.shape[0], mult, -1, p)

    # Quantize
    codes = torch.empty(weight.shape[:-1], device=device, dtype=torch.uint8)
    for i in range(0, weight.shape[0], 16):
        codes[i : i + 16] = torch.argmax(2 * weight[i : i + 16] @ grid.T - grid_norm_2, dim=-1).to(torch.uint8)
    del weight

    codes = codes.reshape(codes.shape[0], -1)
    scales = scales / sqrt(hadamard_size)

    weight, scales, tables, tables2, tune_metadata = prepare_data_transposed(
        codes,
        torch.repeat_interleave(scales.to(dtype), hadamard_size // group_size, dim=1),
        grid.to(dtype),
        num_bits=bits,
        group_size=group_size,
        vector_size=p,
        dtype=dtype,
        device=device,
        check_correctness=False,
    )

    return {
        "weight": weight,
        "scales": scales,
        "tables": tables,
        "tables2": tables2.view(dtype=torch.float16),
        "tune_metadata": tune_metadata,
    }