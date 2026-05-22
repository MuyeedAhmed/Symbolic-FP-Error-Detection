    def test_cutlass_backend_matmul_same_tensor(self):
        max_autotune_gemm_backends = "CUTLASS"

        M = 128
        A = torch.randn(M, M).to(GPU_TYPE).half()

        with config.patch(
            {
                "max_autotune": True,
                "max_autotune_gemm_backends": max_autotune_gemm_backends,
                "cutlass.cutlass_max_profiling_configs": 2,
            }
        ):
            compiled = torch.compile(torch.mm)

            torch.testing.assert_close(A @ A.t(), compiled(A, A.t()))