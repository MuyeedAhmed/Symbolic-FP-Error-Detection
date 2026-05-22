    def test_full_projector_roundtrip_quality(self):
        """project → project_back captures the dominant gradient directions."""
        torch.manual_seed(42)
        # Create a gradient with clear low-rank structure
        u = torch.randn(32, 4)
        v = torch.randn(4, 16)
        grad = u @ v  # rank-4 gradient

        proj = GaLoreProjector(rank = 4, update_proj_gap = 1, scale = 1.0)
        low = proj.project(grad, step = 0)
        reconstructed = proj.project_back(low)

        # For a rank-4 gradient with rank-4 projection, reconstruction
        # should be very close to original
        relative_error = (grad - reconstructed).norm() / grad.norm()
        assert (
            relative_error < 0.05
        ), f"Reconstruction error too high: {relative_error:.4f}"