    def test_inv(self):
        x = np.random.rand(4, 3, 3)
        x_inv = ops.convert_to_numpy(linalg.inv(x))
        x_reconstructed = x @ x_inv
        # high tolerance due to numerical instability
        self.assertAllClose(
            x_reconstructed, np.repeat(np.eye(3)[None], 4, 0), atol=1e-3
        )