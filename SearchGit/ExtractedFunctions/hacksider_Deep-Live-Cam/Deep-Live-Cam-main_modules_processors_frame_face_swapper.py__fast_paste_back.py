def _fast_paste_back(target_img: Frame, bgr_fake: np.ndarray, aimg: np.ndarray, M: np.ndarray) -> Frame:
    """Paste bgr_fake back onto target_img via the inverse affine of M.

    Restricts work to the face bbox in output coordinates and warps a
    precomputed feathered alpha template per-frame instead of running a
    size-scaled erode+blur on the warped mask. Cost is O(crop_area) regardless
    of how much of the frame the face occupies.
    """
    h, w = target_img.shape[:2]
    face_h, face_w = aimg.shape[:2]
    # inswapper's aligned-face space is square (128x128). _get_soft_alpha
    # caches a single NxN template keyed by N, so fail loudly if that ever
    # stops being true rather than silently mis-warping the alpha mask.
    assert face_h == face_w, f"Expected square aligned face, got {face_h}x{face_w}"
    IM = cv2.invertAffineTransform(M)

    # Bbox in output coords from the affine corners of the aligned-face square.
    corners = np.array(
        [[0, 0], [face_w, 0], [face_w, face_h], [0, face_h]], dtype=np.float32
    )
    transformed = (IM[:, :2] @ corners.T).T + IM[:, 2]
    x1 = int(np.floor(transformed[:, 0].min()))
    x2 = int(np.ceil(transformed[:, 0].max()))
    y1 = int(np.floor(transformed[:, 1].min()))
    y2 = int(np.ceil(transformed[:, 1].max()))
    if x1 >= x2 or y1 >= y2:
        return target_img

    # Small interpolation margin only — the feather is baked into the template.
    pad = 2
    y1p, y2p = max(0, y1 - pad), min(h, y2 + pad + 1)
    x1p, x2p = max(0, x1 - pad), min(w, x2 + pad + 1)

    IM_crop = IM.copy()
    IM_crop[0, 2] -= x1p
    IM_crop[1, 2] -= y1p
    crop_w, crop_h = x2p - x1p, y2p - y1p

    soft_alpha = _get_soft_alpha(face_h)
    bgr_fake_crop = cv2.warpAffine(bgr_fake, IM_crop, (crop_w, crop_h), borderMode=cv2.BORDER_REPLICATE)
    alpha_crop = cv2.warpAffine(soft_alpha, IM_crop, (crop_w, crop_h), borderValue=0)

    target_crop = target_img[y1p:y2p, x1p:x2p]

    if _HAS_TORCH_CUDA:
        # Scale alpha to [0, 1] on device — cheaper to upload uint8 than float.
        mask_t = torch.from_numpy(alpha_crop).cuda().float().mul_(1.0 / 255.0).unsqueeze(2)
        fake_t = torch.from_numpy(bgr_fake_crop).float().cuda()
        tgt_t = torch.from_numpy(target_crop).float().cuda()
        blended = (mask_t * fake_t + (1.0 - mask_t) * tgt_t).to(torch.uint8).cpu().numpy()
        target_img[y1p:y2p, x1p:x2p] = blended
    else:
        # Fused uint8 blend via cv2 SIMD — no float32 round-trip.
        # Measured ~7-8× faster than the old numpy float32 path on a 1000×1000 crop.
        alpha_3c = cv2.merge([alpha_crop, alpha_crop, alpha_crop])
        inv_alpha = 255 - alpha_3c
        a_fake = cv2.multiply(bgr_fake_crop, alpha_3c, scale=1.0 / 255.0)
        a_tgt = cv2.multiply(target_crop, inv_alpha, scale=1.0 / 255.0)
        target_img[y1p:y2p, x1p:x2p] = cv2.add(a_fake, a_tgt)

    return target_img