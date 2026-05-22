    def _mkl_matrix(cov_s, cov_r, eps):
        """Compute MKL 3x3 transform matrix from source and ref covariances."""
        eig_val_s, eig_vec_s = torch.linalg.eigh(cov_s)
        sqrt_val_s = torch.sqrt(eig_val_s.clamp_min(0)).clamp_min_(eps)

        scaled_V = eig_vec_s * sqrt_val_s.unsqueeze(0)
        mid = scaled_V.T @ cov_r @ scaled_V
        eig_val_m, eig_vec_m = torch.linalg.eigh(mid)
        sqrt_m = torch.sqrt(eig_val_m.clamp_min(0))

        inv_sqrt_s = 1.0 / sqrt_val_s
        inv_scaled_V = eig_vec_s * inv_sqrt_s.unsqueeze(0)
        M_half = (eig_vec_m * sqrt_m.unsqueeze(0)) @ eig_vec_m.T
        return inv_scaled_V @ M_half @ inv_scaled_V.T