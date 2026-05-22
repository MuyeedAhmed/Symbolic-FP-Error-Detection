def _compute_mask_overlap(masks_a, masks_b):
    """Max of IoU and IoM (intersection over minimum area). More robust to size differences."""
    a_flat = (masks_a > 0).float().flatten(1)
    b_flat = (masks_b > 0).float().flatten(1)
    intersection = a_flat @ b_flat.T
    area_a = a_flat.sum(1, keepdim=True)
    area_b = b_flat.sum(1, keepdim=True).T
    iou = intersection / (area_a + area_b - intersection).clamp(min=1)
    iom = intersection / torch.min(area_a.expand_as(iou), area_b.expand_as(iou)).clamp(min=1)
    return torch.max(iou, iom)