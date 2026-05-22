    def _get_orthogonal_matrices(self, device, dtype):
        """Compute the orthogonal rotation matrices R from BOFT blocks."""
        v = self.weights
        blocks = v[0].to(device=device, dtype=dtype)
        alpha = v[2]
        if alpha is None:
            alpha = 0

        boft_m, block_num, boft_b, _ = blocks.shape
        I = torch.eye(boft_b, device=device, dtype=dtype)

        # Q = blocks - blocks^T (skew-symmetric)
        q = blocks - blocks.transpose(-1, -2)
        normed_q = q

        # Apply constraint if alpha > 0
        if alpha > 0:
            q_norm = torch.norm(q) + 1e-8
            if q_norm > alpha:
                normed_q = q * alpha / q_norm

        # Cayley transform: R = (I + Q)(I - Q)^-1
        r = (I + normed_q) @ (I - normed_q).float().inverse()
        return r, boft_m, boft_b