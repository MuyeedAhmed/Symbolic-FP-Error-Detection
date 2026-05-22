    def test_sp24_matmuls_mat_vec(self) -> None:
        a = torch.randn([64, 128], device="cuda", dtype=torch.float16)
        b = torch.randn([128], device="cuda", dtype=torch.float16)
        a_m = sparse24_largest_mask_2d(a)
        a_s = to_sparse_semi_structured(a)

        with pytest.raises(NotImplementedError):
            torch.testing.assert_close(a_s @ b, (a * a_m) @ b, **atol_rtol_kw[a.dtype])