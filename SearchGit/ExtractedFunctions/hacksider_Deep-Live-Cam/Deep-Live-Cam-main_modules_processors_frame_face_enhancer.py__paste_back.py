def _paste_back(
    frame: Frame,
    enhanced_face: np.ndarray,
    affine_matrix: np.ndarray,
    output_size: int,
) -> Frame:
    """
    Paste an enhanced (aligned) face back onto the original frame using the
    inverse affine transform with feathered-edge blending.

    Optimized: operates on a tight crop around the face bbox instead of the
    full frame, and uses GPU for blending when available.
    """
    h, w = frame.shape[:2]
    inv_matrix = cv2.invertAffineTransform(affine_matrix)

    # Build or reuse cached feathered mask (uint8 — blended via cv2 SIMD ops)
    if _enhancer_cache['mask_size'] != output_size:
        face_mask_f = np.ones((output_size, output_size), dtype=np.float32)
        border = max(1, int(output_size * 0.05))
        ramp_up = np.linspace(0.0, 1.0, border, dtype=np.float32)
        ramp_down = np.linspace(1.0, 0.0, border, dtype=np.float32)
        face_mask_f[:border, :] *= ramp_up[:, None]
        face_mask_f[-border:, :] *= ramp_down[:, None]
        face_mask_f[:, :border] *= ramp_up[None, :]
        face_mask_f[:, -border:] *= ramp_down[None, :]
        _enhancer_cache['mask'] = (face_mask_f * 255.0).astype(np.uint8)
        _enhancer_cache['mask_size'] = output_size

    # Compute tight bbox from affine corners (avoids full-frame warpAffine scan)
    corners = np.array([[0, 0], [output_size, 0],
                        [output_size, output_size], [0, output_size]],
                       dtype=np.float32)
    transformed = (inv_matrix[:, :2] @ corners.T).T + inv_matrix[:, 2]
    x1 = max(0, int(np.floor(transformed[:, 0].min())))
    x2 = min(w, int(np.ceil(transformed[:, 0].max())))
    y1 = max(0, int(np.floor(transformed[:, 1].min())))
    y2 = min(h, int(np.ceil(transformed[:, 1].max())))
    if x1 >= x2 or y1 >= y2:
        return frame

    # Pad a few pixels for feathering
    pad = max(1, int(output_size * 0.05)) + 2
    y1p, y2p = max(0, y1 - pad), min(h, y2 + pad)
    x1p, x2p = max(0, x1 - pad), min(w, x2 + pad)
    crop_w, crop_h = x2p - x1p, y2p - y1p

    # Warp enhanced face and mask into crop space only
    inv_crop = inv_matrix.copy()
    inv_crop[0, 2] -= x1p
    inv_crop[1, 2] -= y1p

    inv_restored_crop = cv2.warpAffine(
        enhanced_face, inv_crop, (crop_w, crop_h),
        borderMode=cv2.BORDER_CONSTANT, borderValue=(0, 0, 0),
    )
    inv_mask_crop = cv2.warpAffine(
        _enhancer_cache['mask'], inv_crop, (crop_w, crop_h),
        borderMode=cv2.BORDER_CONSTANT, borderValue=0,
    )

    target_crop = frame[y1p:y2p, x1p:x2p]

    if _HAS_TORCH_CUDA:
        # Upload uint8 alpha — smaller transfer, scale on device.
        mask_t = torch.from_numpy(inv_mask_crop).cuda().float().mul_(1.0 / 255.0).unsqueeze(2)
        enhanced_t = torch.from_numpy(inv_restored_crop).float().cuda()
        target_t = torch.from_numpy(target_crop).float().cuda()
        blended = (mask_t * enhanced_t + (1.0 - mask_t) * target_t
                   ).to(torch.uint8).cpu().numpy()
        frame[y1p:y2p, x1p:x2p] = blended
    else:
        # Fused uint8 blend via cv2 SIMD — ~7× faster than the float32 round-trip.
        alpha_3c = cv2.merge([inv_mask_crop, inv_mask_crop, inv_mask_crop])
        inv_alpha = 255 - alpha_3c
        a_enh = cv2.multiply(inv_restored_crop, alpha_3c, scale=1.0 / 255.0)
        a_tgt = cv2.multiply(target_crop, inv_alpha, scale=1.0 / 255.0)
        frame[y1p:y2p, x1p:x2p] = cv2.add(a_enh, a_tgt)

    return frame