    def test_sp24_matmuls(self, dtype) -> None:
        M, N, K = 64, 256, 1024
        a = torch.randn([M, K], device="cuda", dtype=dtype)
        b = torch.randn([K, N], device="cuda", dtype=dtype)
        a_m = sparse24_largest_mask_2d(a)
        b_m = sparse24_largest_mask_2d(b)
        (packed, meta, packed_t, meta_t, bitmask) = torch._sparse_semi_structured_tile(
            a
        )
        a_s = SparseSemiStructuredTensorCUTLASS(
            a.shape,
            packed=packed,
            meta=meta,
            packed_t=packed_t,
            meta_t=meta_t,
            compressed_swizzled_bitmask=bitmask,
        )
        (packed, meta, packed_t, meta_t, bitmask) = torch._sparse_semi_structured_tile(
            b
        )
        b_s = SparseSemiStructuredTensorCUTLASS(
            b.shape,
            packed=packed,
            meta=meta,
            packed_t=packed_t,
            meta_t=meta_t,
            compressed_swizzled_bitmask=bitmask,
        )

        torch.testing.assert_close(a_s @ b, (a * a_m) @ b, rtol=1e-1, atol=1.5e-1)
        torch.testing.assert_close(a @ b_s, a @ (b * b_m), rtol=1e-1, atol=1.5e-1)
        torch.testing.assert_close(
            a @ a_s.t(), a @ (a * a_m).t(), rtol=1e-1, atol=1.5e-1
        )
        torch.testing.assert_close(a_s.t() @ a, (a * a_m).t() @ a, rtol=1e-1, atol=1e-1)