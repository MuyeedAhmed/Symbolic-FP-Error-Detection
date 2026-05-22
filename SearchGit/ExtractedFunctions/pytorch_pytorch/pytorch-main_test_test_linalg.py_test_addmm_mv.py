    def test_addmm_mv(self, device, dtype, transpose_a, transpose_b, alpha, beta):
        def gen_mat(w, h, use_transpose: bool = False):
            if not use_transpose:
                return torch.rand(w, h, dtype=dtype, device=device)
            return torch.rand(h, w, dtype=dtype, device=device).t()
        # Regression tests for https://github.com/pytorch/pytorch/issues/136299
        # Should only expose problems on aarch64, but let's be thorough
        m, n , k = 1, 8, 32
        A = gen_mat(m, k, transpose_a)
        B = gen_mat(k, n, transpose_b)
        C = torch.ones(m, n, dtype=dtype, device=device)
        rc = torch.addmm(C, A, B, alpha=alpha, beta=beta)
        ref = alpha * A @ B + beta * C
        self.assertEqual(rc, ref)