        def run_test(shape, batch, nrhs, hermitian):
            A = random_hermitian_pd_matrix(shape, *batch, dtype=dtype, device=device)
            B = make_tensor((*A.shape[:-1], nrhs), dtype=dtype, device=device)
            factors, pivots, info = torch.linalg.ldl_factor_ex(A, hermitian=hermitian)
            X = torch.linalg.ldl_solve(factors, pivots, B, hermitian=hermitian)

            def symmetric(A):
                return A.tril() + A.tril(-1).mT

            # verify A @ X == B
            expected_B = symmetric(A) @ X if not hermitian else A @ X
            self.assertEqual(B, expected_B)