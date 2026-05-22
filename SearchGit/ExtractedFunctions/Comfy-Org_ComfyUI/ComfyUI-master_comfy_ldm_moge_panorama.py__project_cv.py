def _project_cv(points: np.ndarray, extrinsics: np.ndarray, intrinsics: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """World coords -> (uv, depth) in the camera (OpenCV convention)."""
    pts_h = np.concatenate([points, np.ones_like(points[..., :1])], axis=-1)
    cam = pts_h @ extrinsics.T
    cam_xyz = cam[..., :3]
    depth = cam_xyz[..., 2]
    proj = cam_xyz @ intrinsics.T
    uv = proj[..., :2] / proj[..., 2:3].clip(1e-12)
    return uv.astype(np.float32), depth.astype(np.float32)