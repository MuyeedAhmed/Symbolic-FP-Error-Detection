    def execute(cls, splat, translate_x, translate_y, translate_z,
                rotate_x, rotate_y, rotate_z, scale_x, scale_y, scale_z) -> IO.NodeOutput:
        pos = splat.positions
        dev, dt = pos.device, pos.dtype
        q_rot = _euler_to_quat(rotate_x, rotate_y, rotate_z).to(device=dev, dtype=dt)
        R = _quat_to_mat(q_rot[None])[0]                            # (3, 3) node rotation
        D = torch.tensor([scale_x, scale_y, scale_z], dtype=dt, device=dev)
        A = D[:, None] * R                                          # diag(D) @ R: per-axis scale after rotation
        t = torch.tensor([translate_x, translate_y, translate_z], dtype=dt, device=dev)

        positions = pos @ A.T + t                                   # rotate, scale per-axis, then translate
        if scale_x == scale_y == scale_z:                           # uniform: rotation/scale factor out cleanly
            scales = splat.scales * scale_x
            rotations = _quat_mul(q_rot.expand_as(splat.rotations), splat.rotations)
            rotations = rotations / rotations.norm(dim=-1, keepdim=True).clamp_min(1e-12)
        else:                                                       # non-uniform: transform Sigma = A R s^2 R^T A^T, re-extract
            rg = _quat_to_mat(splat.rotations.reshape(-1, 4))       # (M,3,3) per-splat rotation
            s2 = splat.scales.reshape(-1, 3).square()
            cov = (rg * s2[:, None, :]) @ rg.transpose(-1, -2)      # Sigma
            cov = A @ cov @ A.T                                     # A Sigma A^T (A broadcast over splats)
            lam, V = torch.linalg.eigh(cov)                         # symmetric -> eigenvalues (asc), orthonormal axes
            V = V * torch.where(torch.linalg.det(V) < 0, -1.0, 1.0)[..., None, None]   # keep a proper rotation
            scales = lam.clamp_min(0).sqrt().reshape(splat.scales.shape)
            rotations = _mat_to_quat(V).reshape(splat.rotations.shape)
        out = Types.SPLAT(positions, scales, rotations, splat.opacities, splat.sh,
                             counts=getattr(splat, "counts", None))
        return IO.NodeOutput(out)