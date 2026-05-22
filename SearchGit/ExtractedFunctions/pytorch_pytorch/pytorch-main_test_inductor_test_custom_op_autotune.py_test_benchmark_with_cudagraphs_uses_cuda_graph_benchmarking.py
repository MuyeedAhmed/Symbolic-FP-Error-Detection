    def test_benchmark_with_cudagraphs_uses_cuda_graph_benchmarking(self):
        """Test that benchmark_with_cudagraphs flag causes CUDA graph benchmarking to be used."""
        if self.device != "cuda":
            self.skipTest("CUDA graph test requires CUDA device")

        from unittest.mock import patch

        test_op_name = f"test_lib::cudagraph_patch_{id(self)}"

        def fast_decomposition(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight

        @torch.library.custom_op(test_op_name, mutates_args=())
        def cudagraph_patch_op(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight

        @cudagraph_patch_op.register_fake
        def _(x: torch.Tensor, weight: torch.Tensor):
            return torch.empty(
                x.shape[0], weight.shape[1], device=x.device, dtype=x.dtype
            )

        register_custom_op_autotuning(
            cudagraph_patch_op,
            configs=[CustomOpConfig(fast_decomposition)],
            name="cudagraph_patch_autotuned",
            benchmark_with_cudagraphs=True,
            input_gen_fns={
                "x": lambda t: torch.randn_like(t, device=self.device),
                "weight": lambda t: torch.randn_like(t, device=self.device),
            },
        )

        test_x = torch.randn(64, 256, device=self.device, dtype=self.dtype)
        test_weight = torch.randn(256, 128, device=self.device, dtype=self.dtype)

        @torch.compile
        def test_model(x, weight):
            return cudagraph_patch_op(x, weight)

        cuda_graph_benchmark_called = False
        original_benchmark_gpu_with_cuda_graph = torch._inductor.runtime.benchmarking.Benchmarker.benchmark_gpu_with_cuda_graph

        def patched_benchmark_gpu_with_cuda_graph(self, fn):
            nonlocal cuda_graph_benchmark_called
            cuda_graph_benchmark_called = True
            return original_benchmark_gpu_with_cuda_graph(self, fn)

        torch._dynamo.reset()
        with config.patch(max_autotune=True, fx_graph_cache=False):
            with patch.object(
                torch._inductor.runtime.benchmarking.Benchmarker,
                "benchmark_gpu_with_cuda_graph",
                patched_benchmark_gpu_with_cuda_graph,
            ):
                result = test_model(test_x, test_weight)

        self.assertTrue(
            cuda_graph_benchmark_called,
            "benchmark_gpu_with_cuda_graph should have been called",
        )
        torch.testing.assert_close(result, test_x @ test_weight, rtol=1e-1, atol=1e-1)