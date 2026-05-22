def _unproject_cv(uv: np.ndarray, depth: np.ndarray,
                  extrinsics: np.ndarray, intrinsics: np.ndarray) -> np.ndarray:
    """Back-project pixels into world coords (OpenCV convention)."""
    pix = np.concatenate([uv, np.ones_like(uv[..., :1])], axis=-1)
    K_inv = np.linalg.inv(intrinsics)
    cam = pix @ K_inv.T * depth[..., None]
    cam_h = np.concatenate([cam, np.ones_like(cam[..., :1])], axis=-1)
    E_inv = np.linalg.inv(extrinsics)
    return (cam_h @ E_inv.T)[..., :3]