def depth_map_to_point_map(depth: torch.Tensor, intrinsics: torch.Tensor) -> torch.Tensor:
    """Back-project a (..., H, W) depth map through K^-1 to (..., H, W, 3) camera-space points.

    Intrinsics use normalized image coords (x in [0, 1] left->right, y in [0, 1] top->bottom).
    """
    H, W = depth.shape[-2:]
    device, dtype = depth.device, depth.dtype
    u = (torch.arange(W, dtype=dtype, device=device) + 0.5) / W
    v = (torch.arange(H, dtype=dtype, device=device) + 0.5) / H
    grid_v, grid_u = torch.meshgrid(v, u, indexing="ij")
    pix = torch.stack([grid_u, grid_v, torch.ones_like(grid_u)], dim=-1)
    K_inv = torch.linalg.inv(intrinsics)
    rays = torch.einsum("...ij,hwj->...hwi", K_inv, pix)
    return rays * depth.unsqueeze(-1)