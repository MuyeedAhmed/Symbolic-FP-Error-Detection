    def test_partial_input_gen_fns(self):
        """Test autotuning when input_gen_fns covers only some inputs.

        The uncovered inputs should fall back to ir_node_to_tensor with concrete hints.
        """
        test_op_name = f"test_lib::partial_gen_{id(self)}"

        def decomposition(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight

        @torch.library.custom_op(test_op_name, mutates_args=())
        def partial_gen_op(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight

        @partial_gen_op.register_fake
        def _(x: torch.Tensor, weight: torch.Tensor):
            return torch.empty(
                x.shape[0], weight.shape[1], device=x.device, dtype=x.dtype
            )

        # Only provide input_gen_fn for "x", not "weight"
        register_custom_op_autotuning(
            partial_gen_op,
            configs=[CustomOpConfig(decomposition)],
            name="partial_gen_autotuned",
            input_gen_fns={
                "x": lambda t: torch.randn_like(t, device=self.device),
            },
        )

        test_x = torch.randn(64, 256, device=self.device, dtype=self.dtype)
        test_weight = torch.randn(256, 128, device=self.device, dtype=self.dtype)

        @torch.compile
        def test_model(x, weight):
            return partial_gen_op(x, weight)

        torch._dynamo.reset()
        with config.patch(max_autotune=True, fx_graph_cache=False):
            result = test_model(test_x, test_weight)

        torch.testing.assert_close(result, test_x @ test_weight, rtol=1e-1, atol=1e-1)