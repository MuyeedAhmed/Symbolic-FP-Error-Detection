    def test_min_speedup_threshold_api(self):
        """Test that min_speedup_threshold parameter is accepted and compilation works."""
        test_op_name = f"test_lib::min_speedup_{id(self)}"

        def decomposition(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight

        @torch.library.custom_op(test_op_name, mutates_args=())
        def min_speedup_op(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight

        @min_speedup_op.register_fake
        def _(x: torch.Tensor, weight: torch.Tensor):
            return torch.empty(
                x.shape[0], weight.shape[1], device=x.device, dtype=x.dtype
            )

        # Test that API accepts min_speedup_threshold parameter
        register_custom_op_autotuning(
            min_speedup_op,
            configs=[CustomOpConfig(decomposition)],
            name="min_speedup_autotuned",
            min_speedup_threshold=1.5,
            input_gen_fns={
                "x": lambda t: torch.randn_like(t, device=self.device),
                "weight": lambda t: torch.randn_like(t, device=self.device),
            },
        )

        test_x = torch.randn(64, 256, device=self.device, dtype=self.dtype)
        test_weight = torch.randn(256, 128, device=self.device, dtype=self.dtype)

        @torch.compile
        def test_model(x, weight):
            return min_speedup_op(x, weight)

        torch._dynamo.reset()
        with config.patch(max_autotune=True, fx_graph_cache=False):
            result = test_model(test_x, test_weight)

        torch.testing.assert_close(result, test_x @ test_weight, rtol=1e-1, atol=1e-1)