    def test_blockwise_nvfp4_compile(self, device) -> None:

        M, K, N = 128, 128, 128
        BLOCK_SIZE = 32 if torch.version.hip else 16
        fp4_scaling_dtype = torch.float8_e8m0fnu if torch.version.hip else torch.float8_e4m3fn

        A_ref = torch.eye(M, device=device, dtype=torch.bfloat16)
        B_ref = torch.eye(M, device=device, dtype=torch.bfloat16)

        A = _bfloat16_to_float4_e2m1fn_x2(A_ref)
        B = _bfloat16_to_float4_e2m1fn_x2(B_ref)

        A_scale = torch.full((M, ceil_div(K, BLOCK_SIZE)), 1.0, device=device, dtype=fp4_scaling_dtype)
        B_scale = torch.full((N, ceil_div(K, BLOCK_SIZE)), 1.0, device=device, dtype=fp4_scaling_dtype)
        C_ref = A_ref @ B_ref.t()

        compiled_scaled_mm = torch.compile(scaled_mm_wrap, backend="inductor")
        # C = scaled_mm_wrap(
        C = compiled_scaled_mm(
            A,
            B.t(),
            A_scale,
            B_scale,
            out_dtype=torch.bfloat16,
            use_fast_accum=False,
        )
        torch.testing.assert_close(C, C_ref, atol=0, rtol=0)