    def _test_linalg_solve_triangular(self, A, B, upper, left, uni):
        X = torch.linalg.solve_triangular(A, B, upper=upper, left=left, unitriangular=uni)
        if left:
            self.assertEqual(A @ X, B)
        else:
            self.assertEqual(X @ A, B)
        out = B
        # B may be expanded
        if not B.is_contiguous() and not B.transpose(-2, -1).is_contiguous():
            out = B.clone()
        torch.linalg.solve_triangular(A, B, upper=upper, left=left, unitriangular=uni, out=out)
        self.assertEqual(X, out)