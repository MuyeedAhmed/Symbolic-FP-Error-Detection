    def test_cutlass_backend_matmul_nonzero_offset(self):
        max_autotune_gemm_backends = "CUTLASS"

        M = 129
        A = torch.randn(M, M - 1).to(GPU_TYPE).half()

        with config.patch(
            {
                "max_autotune": True,
                "max_autotune_gemm_backends": max_autotune_gemm_backends,
                "cutlass.cutlass_max_profiling_configs": 2,
            }
        ):
            compiled = torch.compile(torch.mm)
            torch.testing.assert_close(
                A[1:, :] @ A[1:, :].t(), compiled(A[1:, :], A[1:, :].t())
            )