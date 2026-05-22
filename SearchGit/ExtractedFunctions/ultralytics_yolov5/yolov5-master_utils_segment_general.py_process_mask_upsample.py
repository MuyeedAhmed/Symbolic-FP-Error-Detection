def process_mask_upsample(protos, masks_in, bboxes, shape):
    """Crop after upsample.

    Args:
        protos: [mask_dim, mask_h, mask_w]
        masks_in: [n, mask_dim], n is number of masks after nms
        bboxes: [n, 4], n is number of masks after nms
        shape: input_image_size, (h, w).

    Returns:
        h, w, n
    """
    c, mh, mw = protos.shape  # CHW
    masks = (masks_in @ protos.float().view(c, -1)).sigmoid().view(-1, mh, mw)
    masks = F.interpolate(masks[None], shape, mode="bilinear", align_corners=False)[0]  # CHW
    masks = crop_mask(masks, bboxes)  # CHW
    return masks.gt_(0.5)