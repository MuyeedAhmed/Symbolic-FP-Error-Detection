    def test_linalg_solve_triangular_broadcasting(self, device, dtype):
        make_arg = partial(make_tensor, dtype=dtype, device=device)

        sizes = (((2, 1, 3, 4, 4), (2, 1, 3, 4, 6)),
                 ((2, 1, 3, 4, 4), (4, 6)),
                 ((4, 4), (2, 1, 3, 4, 2)),
                 ((1, 3, 1, 4, 4), (2, 1, 3, 4, 5)))
        for size_A, size_B in sizes:
            for left, upper, uni in itertools.product([True, False], repeat=3):
                A = make_arg(size_A)
                if upper:
                    A.triu_()
                else:
                    A.tril_()
                diag = A.diagonal(0, -2, -1)
                if uni:
                    diag.fill_(1.)
                else:
                    diag[diag.abs() < 1e-6] = 1.
                B = make_arg(size_B)
                if not left:
                    B.transpose_(-2, -1)

                X = torch.linalg.solve_triangular(A, B, upper=upper, left=left, unitriangular=uni)
                if left:
                    B_other = A @ X
                else:
                    B_other = X @ A

                self.assertEqual(*torch.broadcast_tensors(B, B_other))