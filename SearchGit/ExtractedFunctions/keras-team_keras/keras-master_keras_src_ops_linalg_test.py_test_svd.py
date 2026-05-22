    def test_svd(self):
        x = np.random.rand(4, 30, 20).astype("float32")
        u, s, vh = linalg.svd(x)
        x_reconstructed = (u[..., :, : s.shape[-1]] * s[..., None, :]) @ vh[
            ..., : s.shape[-1], :
        ]
        # High tolerance due to numerical instability
        self.assertAllClose(
            x_reconstructed, x, atol=1e-3, tpu_atol=1e-2, tpu_rtol=1e-2
        )

        # Test `compute_uv=False`
        s_no_uv = linalg.svd(x, compute_uv=False)
        self.assertAllClose(
            s_no_uv, s, atol=1e-5, rtol=1e-5, tpu_atol=1e-2, tpu_rtol=1e-2
        )