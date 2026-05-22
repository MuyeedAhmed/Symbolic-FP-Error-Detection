    def test_fp16_backward_preserves_subnormals_rocm(self, backend, batched):
        # Regression test for issue #182952. On ROCm, the hipBLASLt path for
        # at::Half had no equivalent of rocBLAS's fp16_alt_impl, so backward
        # fp16 GEMMs silently flushed subnormals to zero. The dispatcher now
        # routes fp16 backward GEMMs to rocBLAS on gfx90a regardless of the
        # user's preferred backend.
        dtype = torch.float16
        M = K = N = 64
        sub_val = torch.finfo(dtype).tiny / 2
        ref = M * sub_val
        with blas_library_context(backend):
            if batched:
                B = 3
                x = torch.ones(B, M, K, dtype=dtype, device="cuda")
                w = torch.nn.Parameter(
                    torch.ones(B, K, N, dtype=dtype, device="cuda")
                )
                d = torch.full((B, M, N), sub_val, dtype=dtype, device="cuda")
            else:
                x = torch.ones(M, K, dtype=dtype, device="cuda")
                w = torch.nn.Parameter(
                    torch.ones(K, N, dtype=dtype, device="cuda")
                )
                d = torch.full((M, N), sub_val, dtype=dtype, device="cuda")
            (x @ w).backward(d)
            torch.cuda.synchronize()
            grad = w.grad.flatten()[0].item()
            self.assertEqual(grad, ref)