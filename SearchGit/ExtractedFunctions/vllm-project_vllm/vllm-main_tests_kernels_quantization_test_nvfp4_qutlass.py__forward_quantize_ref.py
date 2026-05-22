def _forward_quantize_ref(x: torch.Tensor, h: torch.Tensor, rot_size: int):
    device = x.device

    xh_ref64 = (
        x.unflatten(dim=-1, sizes=(-1, rot_size)).to(dtype=torch.float64)
        @ h.reshape(rot_size, rot_size).to(dtype=torch.float64)
    ).flatten(start_dim=-2)

    abs_max = xh_ref64.unflatten(dim=-1, sizes=(-1, 16)).abs().amax(dim=-1)
    scales_ref64_ = abs_max + 1e-8

    xh_e4m3_ref = scales_ref64_.to(dtype=torch.float8_e4m3fn)
    scales_ref64 = xh_e4m3_ref.to(dtype=torch.float64)
    xh_scaled_ref64 = (
        xh_ref64.unflatten(dim=-1, sizes=(-1, 16)) / scales_ref64[..., None]
    ).flatten(start_dim=-2)

    xh_scaled_ref64 *= 6.0

    clip_mask_unpacked_ref = xh_scaled_ref64.abs() < 6.0
    clip_mask_ref = torch.zeros(
        *x.shape[:-1], x.size(-1) // 8, dtype=torch.uint8, device=device
    )
    for i in range(8):
        clip_mask_ref |= clip_mask_unpacked_ref[..., i::8].to(dtype=torch.uint8) << i

    xh_fp4_ref, xh_e2m1_ref = _rtne_fp4(xh_scaled_ref64)
    xh_dq, xh_fp4_dq, scales_dq = _dq_fp4(xh_e2m1_ref, xh_e4m3_ref, 6.0)
    clip_mask_unpacked_dq = _unpack_mask(clip_mask_ref)

    assert xh_fp4_dq.equal(xh_fp4_ref)
    assert scales_dq.equal(scales_ref64)
    assert clip_mask_unpacked_dq.equal(clip_mask_unpacked_ref)

    return (
        xh_dq,
        clip_mask_unpacked_ref,
        (xh_e2m1_ref, xh_e4m3_ref, clip_mask_ref),
    )