    def test_vmap_multi_dot_failure_1D_input(self):
        # special exception for first and last tensors so making giving 3 items avoids special cases
        inputs = (torch.randn(3, 3), torch.randn(3), torch.randn(3, 3))
        with self.assertRaisesRegex(RuntimeError, "tensor 1 must be 2D but got 1D"):
            torch.linalg.multi_dot(inputs)

        # square inputs are more likely to pass linalg checks
        inputs = tuple(i.expand(i.shape[0], i.shape[0]) for i in inputs)
        with self.assertRaisesRegex(RuntimeError, "tensor 1 must be 2D but got 1D"):
            return vmap(torch.linalg.multi_dot)(inputs)