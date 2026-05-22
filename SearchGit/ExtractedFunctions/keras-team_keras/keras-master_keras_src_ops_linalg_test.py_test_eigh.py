    def test_eigh(self):
        x = np.random.rand(2, 3, 3)
        x = x @ x.transpose((0, 2, 1))
        w, v = map(ops.convert_to_numpy, linalg.eigh(x))
        x_reconstructed = (v * w[..., None, :]) @ v.transpose((0, 2, 1))
        self.assertAllClose(x_reconstructed, x, atol=1e-4)