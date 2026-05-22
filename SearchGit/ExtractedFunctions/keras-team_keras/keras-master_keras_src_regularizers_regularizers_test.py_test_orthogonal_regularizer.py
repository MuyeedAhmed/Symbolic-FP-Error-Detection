    def test_orthogonal_regularizer(self):
        value = np.random.random((4, 4)).astype(np.float32)
        x = backend.Variable(value)
        y = regularizers.OrthogonalRegularizer(factor=0.1, mode="rows")(x)

        l2_norm = np.linalg.norm(value, axis=1, keepdims=True)
        inputs = value / l2_norm
        self.assertAllClose(
            y,
            0.1
            * 0.5
            * np.sum(
                np.abs(np.dot(inputs, np.transpose(inputs)) * (1.0 - np.eye(4)))
            )
            / (4.0 * (4.0 - 1.0) / 2.0),
            tpu_atol=1e-4,
            tpu_rtol=1e-4,
        )