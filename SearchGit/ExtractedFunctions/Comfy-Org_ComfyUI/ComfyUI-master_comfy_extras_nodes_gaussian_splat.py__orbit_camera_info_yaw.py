def _orbit_camera_info_yaw(camera_info, angle_deg, dev):
    # Turntable: rigidly rotate a camera_info about world +Y around its target by angle_deg. Returns a new dict.
    a = math.radians(angle_deg)
    ca, sa = math.cos(a), math.sin(a)
    v = lambda d: torch.tensor([float(d.get("x", 0.0)), float(d.get("y", 0.0)), float(d.get("z", 0.0))], device=dev)
    pos, tgt = v(camera_info.get("position", {})), v(camera_info.get("target", {}))
    Ry = torch.tensor([[ca, 0.0, sa], [0.0, 1.0, 0.0], [-sa, 0.0, ca]], device=dev)
    new_pos = tgt + Ry @ (pos - tgt)
    q = camera_info.get("quaternion") or {}
    qcur = torch.tensor([float(q.get("w", 1.0)), float(q.get("x", 0.0)),
                         float(q.get("y", 0.0)), float(q.get("z", 0.0))], device=dev)
    qy = torch.tensor([math.cos(a / 2), 0.0, math.sin(a / 2), 0.0], device=dev)   # world +Y rotation
    qn = _quat_mul(qy[None], qcur[None])[0]
    xyz = lambda t: {"x": float(t[0]), "y": float(t[1]), "z": float(t[2])}
    return {**camera_info, "position": xyz(new_pos),
            "quaternion": {"x": float(qn[1]), "y": float(qn[2]), "z": float(qn[3]), "w": float(qn[0])}}