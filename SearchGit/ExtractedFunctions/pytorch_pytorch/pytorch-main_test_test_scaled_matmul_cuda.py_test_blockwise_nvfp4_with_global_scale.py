    def test_blockwise_nvfp4_with_global_scale(self, mkn) -> None:
        device = 'cuda'
        M, K, N = mkn
        BLOCK_SIZE = 16
        # Note: SQNR target from `test_blockwise_mxfp8_nvfp4_mxfp4_numerics` test
        approx_match_sqnr_target = 15.8

        A_ref = torch.randn((M, K), device=device, dtype=torch.bfloat16) * 1000
        B_ref = torch.randn((N, K), device=device, dtype=torch.bfloat16) * 1000

        A, A_scale, A_global_scale = data_to_nvfp4_with_global_scale(A_ref, BLOCK_SIZE)
        B, B_scale, B_global_scale = data_to_nvfp4_with_global_scale(B_ref, BLOCK_SIZE)

        if torch.version.cuda:
            A_scale = to_blocked(A_scale)
            B_scale = to_blocked(B_scale)
            swizzle = [SwizzleType.SWIZZLE_32_4_4, SwizzleType.NO_SWIZZLE]
        else:
            swizzle = [SwizzleType.NO_SWIZZLE, SwizzleType.NO_SWIZZLE]

        C_ref = A_ref @ B_ref.t()

        C = scaled_mm(
            A,
            B.t(),
            scale_a=[A_scale, A_global_scale],
            scale_recipe_a=[ScalingType.BlockWise1x16, ScalingType.TensorWise],
            scale_b=[B_scale, B_global_scale],
            scale_recipe_b=[ScalingType.BlockWise1x16, ScalingType.TensorWise],
            swizzle_a=swizzle,
            swizzle_b=swizzle,
            output_dtype=torch.bfloat16,
        )

        sqnr = compute_error(C_ref, C)
        if sqnr.item() <= approx_match_sqnr_target:
            raise AssertionError(f"sqnr {sqnr.item()} should be > {approx_match_sqnr_target}")