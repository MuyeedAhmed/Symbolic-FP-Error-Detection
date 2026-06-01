def _lookat_quat_wxyz(position, target, dev):
    # three.js lookAt in world frame: camera local +Z = (eye - target), up = world +Y. Returns wxyz.
    z = position - target
    z = z / z.norm().clamp_min(1e-8)
    up0 = torch.tensor([0.0, 1.0, 0.0], device=dev)
    if z.dot(up0).abs() > 0.999:                                # looking straight up/down
        up0 = torch.tensor([0.0, 0.0, 1.0], device=dev)
    x = torch.linalg.cross(up0, z)
    x = x / x.norm().clamp_min(1e-8)
    y = torch.linalg.cross(z, x)
    R = torch.stack([x, y, z], dim=1)                           # columns = camera world axes
    return _mat_to_quat(R[None])[0]