    def test_blockwise_mxfp8_compile(self, device) -> None:

        M, K, N = 128, 128, 128
        BLOCK_SIZE = 32

        A_ref = torch.eye(M, device=device, dtype=torch.bfloat16)
        B_ref = torch.eye(M, device=device, dtype=torch.bfloat16)

        A = A_ref.to(torch.float8_e4m3fn)
        B = B_ref.to(torch.float8_e4m3fn)

        A_scale = torch.full((M, ceil_div(K, BLOCK_SIZE)), 1.0, device=device, dtype=torch.float8_e8m0fnu)
        B_scale = torch.full((N, ceil_div(K, BLOCK_SIZE)), 1.0, device=device, dtype=torch.float8_e8m0fnu)
        C_ref = A_ref @ B_ref.t()

        compiled_scaled_mm = torch.compile(scaled_mm_wrap, backend="inductor")
        C = compiled_scaled_mm(
            A,
            B.t(),
            A_scale,
            B_scale,
            out_dtype=torch.bfloat16,
            use_fast_accum=False,
        )
        torch.testing.assert_close(C, C_ref, atol=0, rtol=0)