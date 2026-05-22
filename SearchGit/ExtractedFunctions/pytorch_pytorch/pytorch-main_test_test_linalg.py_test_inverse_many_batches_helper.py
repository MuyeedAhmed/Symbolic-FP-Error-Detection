        def test_inverse_many_batches_helper(torch_inverse, b, n):
            matrices = make_arg(b, n, n)
            matrices_inverse = torch_inverse(matrices)

            # Compare against NumPy output
            expected = np.linalg.inv(matrices.cpu().numpy())
            self.assertEqual(matrices_inverse, expected, atol=self.precision, rtol=1e-3)