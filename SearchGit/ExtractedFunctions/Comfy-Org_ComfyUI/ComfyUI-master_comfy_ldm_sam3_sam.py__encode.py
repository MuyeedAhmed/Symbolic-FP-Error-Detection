    def _encode(self, normalized_coords):
        """Map normalized [0,1] coordinates to fourier features via random projection. Computes in fp32."""
        orig_dtype = normalized_coords.dtype
        proj_matrix = self.positional_encoding_gaussian_matrix.to(device=normalized_coords.device, dtype=torch.float32)
        projected = 2 * math.pi * (2 * normalized_coords.float() - 1) @ proj_matrix
        return torch.cat([projected.sin(), projected.cos()], dim=-1).to(orig_dtype)