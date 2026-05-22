    def test_cublas_batch_invariance_blackwell(self, device, dtype):
        orig_bf16 = torch.backends.cuda.matmul.allow_bf16_reduced_precision_reduction
        orig_fp16 = torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction
        torch.backends.cuda.matmul.allow_bf16_reduced_precision_reduction = (False, False)
        torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction = (False, False)
        with blas_library_context('cublaslt'):
            N = 2048
            K = 6144
            M_max = 32
            x = torch.randn(M_max, K, device="cuda", dtype=torch.bfloat16)
            w = torch.randn(N, K, device="cuda", dtype=torch.bfloat16).t()
            full = x @ w
            xx = x[:1]
            out = xx @ w
            self.assertEqual(full[:1], out, atol=0., rtol=0.)
        torch.backends.cuda.matmul.allow_bf16_reduced_precision_reduction = orig_bf16
        torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction = orig_fp16