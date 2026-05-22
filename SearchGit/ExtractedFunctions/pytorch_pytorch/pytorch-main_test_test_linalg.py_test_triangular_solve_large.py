    def test_triangular_solve_large(self, device, dtype):
        # Repro for https://github.com/pytorch/pytorch/issues/79191
        A = torch.randn(1, 2, 2, device=device, dtype=dtype).tril_()
        B = torch.randn(1, 2, 524281, device=device, dtype=dtype)
        X = torch.linalg.solve_triangular(A, B, upper=False)
        self.assertEqual(A @ X, B)