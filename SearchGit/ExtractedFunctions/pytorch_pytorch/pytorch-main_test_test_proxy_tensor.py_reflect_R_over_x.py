        def reflect_R_over_x(R):
            reflect = torch.eye(3, device=R.device)
            reflect[0, 0] = -1
            return reflect @ R @ reflect