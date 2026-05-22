    def test_pack_both_ways_meta_correctness(self, dtype, backend) -> None:
        M, N = 128, 256
        # Construct x to make sure we always have exactly 8 elements per 4x4 tile
        a = (4 * torch.arange(8))[:, None] + torch.arange(8)[None, :]
        a = a.repeat(M // 8, N // 8)
        if a.shape != (M, N):
            raise AssertionError(f"a.shape should be ({M}, {N}), got {a.shape}")
        a = a.cuda().to(dtype)
        b = torch.randn([a.shape[1], 128], device="cuda", dtype=dtype)

        a_sparse = SEMI_STRUCTURED_SUPPORTED_BACKENDS[backend].prune_dense_static_sort(
            a
        )

        mask_dense = sparse24_largest_mask_2d(a).to(dtype)

        if backend == "cutlass":
            if not isinstance(a_sparse, SparseSemiStructuredTensorCUTLASS):
                raise AssertionError(
                    f"a_sparse should be SparseSemiStructuredTensorCUTLASS, got {type(a_sparse)}"
                )
            (packed, meta, packed_t, meta_t, bitmask) = (
                torch._sparse_semi_structured_tile(mask_dense, use_cutlass=True)
            )

            sparse_mask = SparseSemiStructuredTensorCUTLASS(
                mask_dense.shape,
                packed=packed,
                meta=meta,
                packed_t=packed_t,
                meta_t=meta_t,
                compressed_swizzled_bitmask=bitmask,
            )
            torch.testing.assert_close(
                a_sparse.meta.view(torch.short), sparse_mask.meta
            )

        ref_gemm = (mask_dense * a) @ b
        pack_gemm = a_sparse @ b
        torch.testing.assert_close(ref_gemm, pack_gemm, **atol_rtol_kw[dtype])