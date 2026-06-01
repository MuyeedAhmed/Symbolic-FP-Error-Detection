def _view_matrix_t(yaw_deg, pitch_deg, device):
    y, p = math.radians(yaw_deg), math.radians(pitch_deg)
    cy, sy, cp, sp = math.cos(y), math.sin(y), math.cos(p), math.sin(p)
    Ry = torch.tensor([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]], device=device)
    Rx = torch.tensor([[1, 0, 0], [0, cp, -sp], [0, sp, cp]], device=device)
    return Rx @ Ry