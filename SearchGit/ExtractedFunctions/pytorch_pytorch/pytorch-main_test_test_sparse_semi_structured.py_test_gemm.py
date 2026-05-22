    def test_gemm(self, dtype) -> None:
        M, N, K = 32, 32, 64
        a = torch.randn([M, K], device="cuda", dtype=dtype)
        b = torch.randn([K, N], device="cuda", dtype=dtype)
        mask = rand_sparse_semi_structured_mask(M, K, dtype=torch.bool)

        a.masked_fill_(~mask, 0)

        a_sparse = to_sparse_semi_structured(a)

        masked_a = a * mask
        ref_out = masked_a @ b
        sp24_out = a_sparse @ b
        torch.testing.assert_close(ref_out, sp24_out, **atol_rtol_kw[dtype])