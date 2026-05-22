    def test_hadamard_orthonormal(self, dim):
        """H must be orthonormal: H @ H^T = I."""
        H = _build_hadamard(dim, DEVICE_TYPE)
        eye = H @ H.T
        assert torch.allclose(eye, torch.eye(dim, device=DEVICE_TYPE), atol=1e-5), (
            f"Hadamard not orthonormal for dim={dim}"
        )