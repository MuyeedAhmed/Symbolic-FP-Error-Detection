    def test_orthogonal(self):
        shape = (5, 5)
        gain = 2.0
        seed = 1234
        initializer = initializers.Orthogonal(gain=gain, seed=seed)
        values = initializer(shape=shape)
        self.assertEqual(initializer.seed, seed)
        self.assertEqual(initializer.gain, gain)

        self.assertEqual(values.shape, shape)
        array = backend.convert_to_numpy(values)
        # Making sure that the columns have gain * unit norm value
        for column in array.T:
            self.assertAlmostEqual(np.linalg.norm(column), gain * 1.0)

        # Making sure that each column is orthonormal to the other column
        for i in range(array.shape[-1]):
            for j in range(i + 1, array.shape[-1]):
                self.assertAlmostEqual(
                    np.dot(array[..., i], array[..., j]), 0.0
                )

        self.run_class_serialization_test(initializer)

        # Test compatible class_name
        initializer = initializers.get("OrthogonalInitializer")
        self.assertIsInstance(initializer, initializers.Orthogonal)