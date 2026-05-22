def post_dark_unbiased_data_processing(coords: np.ndarray, batch_heatmaps: np.ndarray, kernel: int = 3) -> np.ndarray:
    """DARK post-pocessing. Implemented by unbiased_data_processing.

    Paper references:
    - Huang et al. The Devil is in the Details: Delving into Unbiased Data Processing for Human Pose Estimation (CVPR 2020).
    - Zhang et al. Distribution-Aware Coordinate Representation for Human Pose Estimation (CVPR 2020).

    Args:
        coords (`np.ndarray` of shape `(num_persons, num_keypoints, 2)`):
            Initial coordinates of human pose.
        batch_heatmaps (`np.ndarray` of shape `(batch_size, num_keypoints, height, width)`):
            Batched heatmaps as predicted by the model.
            A batch_size of 1 is used for the bottom up paradigm where all persons share the same heatmap.
            A batch_size of `num_persons` is used for the top down paradigm where each person has its own heatmaps.
        kernel (`int`, *optional*, defaults to 3):
            Gaussian kernel size (K) for modulation.

    Returns:
        `np.ndarray` of shape `(num_persons, num_keypoints, 2)` ):
            Refined coordinates.
    """
    batch_size, num_keypoints, height, width = batch_heatmaps.shape
    num_coords = coords.shape[0]
    if not (batch_size == 1 or batch_size == num_coords):
        raise ValueError("The batch size of heatmaps should be 1 or equal to the batch size of coordinates.")
    radius = int((kernel - 1) // 2)
    batch_heatmaps = np.array(
        [
            [gaussian_filter(heatmap, sigma=0.8, radius=(radius, radius), axes=(0, 1)) for heatmap in heatmaps]
            for heatmaps in batch_heatmaps
        ]
    )
    batch_heatmaps = np.clip(batch_heatmaps, 0.001, 50)
    batch_heatmaps = np.log(batch_heatmaps)

    batch_heatmaps_pad = np.pad(batch_heatmaps, ((0, 0), (0, 0), (1, 1), (1, 1)), mode="edge").flatten()
    index = coords[..., 0] + 1 + (coords[..., 1] + 1) * (width + 2)
    index += (width + 2) * (height + 2) * np.arange(0, batch_size * num_keypoints).reshape(-1, num_keypoints)
    index = index.astype(int).reshape(-1, 1)
    i_ = batch_heatmaps_pad[index]
    ix1 = batch_heatmaps_pad[index + 1]
    iy1 = batch_heatmaps_pad[index + width + 2]
    ix1y1 = batch_heatmaps_pad[index + width + 3]
    ix1_y1_ = batch_heatmaps_pad[index - width - 3]
    ix1_ = batch_heatmaps_pad[index - 1]
    iy1_ = batch_heatmaps_pad[index - 2 - width]
    dx = 0.5 * (ix1 - ix1_)
    dy = 0.5 * (iy1 - iy1_)
    derivative = np.concatenate([dx, dy], axis=1)
    derivative = derivative.reshape(num_coords, num_keypoints, 2, 1)
    dxx = ix1 - 2 * i_ + ix1_
    dyy = iy1 - 2 * i_ + iy1_
    dxy = 0.5 * (ix1y1 - ix1 - iy1 + i_ + i_ - ix1_ - iy1_ + ix1_y1_)
    hessian = np.concatenate([dxx, dxy, dxy, dyy], axis=1)
    hessian = hessian.reshape(num_coords, num_keypoints, 2, 2)
    hessian = np.linalg.inv(hessian + np.finfo(np.float32).eps * np.eye(2))
    coords -= np.einsum("ijmn,ijnk->ijmk", hessian, derivative).squeeze()
    return coords