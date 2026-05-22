    def test_old_cholesky_batched_upper(self, device, dtype):
        from torch.testing._internal.common_utils import random_hermitian_pd_matrix

        batchsize = 2
        A = random_hermitian_pd_matrix(3, batchsize, dtype=dtype, device=device)
        A_triu = A.triu()  # fill the lower triangular part with zero

        U = torch.cholesky(A_triu, upper=True)

        reconstruct_A = U.mH @ U
        self.assertEqual(A, reconstruct_A)