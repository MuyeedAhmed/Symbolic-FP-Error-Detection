    def test_cholesky_upper_reconstructs(self, device, dtype):
        batch_dims = (1,)
        matrix_size = 65
        A = torch.randn(
            *(batch_dims + (matrix_size, matrix_size)), dtype=dtype, device=device
        )
        pd_matrix = A @ A.mH + torch.eye(matrix_size, dtype=dtype, device=device)
        pd_matrix = pd_matrix.squeeze(0)
        U = torch.linalg.cholesky(pd_matrix, upper=True)
        self.assertEqual(U, torch.triu(U))
        reconstructed = U.mH @ U
        self.assertEqual(pd_matrix, reconstructed, atol=1e-4, rtol=1e-5)