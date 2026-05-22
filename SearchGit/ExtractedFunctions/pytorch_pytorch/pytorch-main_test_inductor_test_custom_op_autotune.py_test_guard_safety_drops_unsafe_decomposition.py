    def test_guard_safety_drops_unsafe_decomposition(self):
        """Test that decompositions adding guards are replaced with fallback.

        Compiles with m=8 (divisible by 4, satisfies the unsafe impl's guard),
        then calls with m=7 (not divisible by 4). If the Mod guard leaked,
        the m=7 call would either crash or produce incorrect results because
        the compiled code assumes m is divisible by 4.
        """
        test_op_name = f"test_lib::guard_safety_{id(self)}"

        def unsafe_impl(mat1: torch.Tensor, mat2: torch.Tensor) -> torch.Tensor:
            m = mat1.shape[0]
            if m % 4 == 0:
                n = mat2.shape[1]
                k = mat1.shape[1]
                m_parts = m // 4
                a = mat1.reshape(4, m_parts, k)
                result = torch.bmm(a, mat2.unsqueeze(0).expand(4, -1, -1))
                return result.reshape(m, n)
            return mat1 @ mat2

        def safe_impl(mat1: torch.Tensor, mat2: torch.Tensor) -> torch.Tensor:
            return mat1 @ mat2

        @torch.library.custom_op(test_op_name, mutates_args=())
        def guard_safety_op(mat1: torch.Tensor, mat2: torch.Tensor) -> torch.Tensor:
            return mat1 @ mat2

        @guard_safety_op.register_fake
        def _(mat1: torch.Tensor, mat2: torch.Tensor):
            return torch.empty(
                mat1.shape[0], mat2.shape[1], device=mat1.device, dtype=mat1.dtype
            )

        register_custom_op_autotuning(
            guard_safety_op,
            configs=[CustomOpConfig(unsafe_impl)],
            name="guard_safety_autotuned",
            dispatch_on={"tensor_name": "mat1", "dim": 0},
            split_points=[4, 16],
            input_gen_fns={
                "mat1": lambda t: torch.randn_like(t, device=self.device) * 0.1,
                "mat2": lambda t: torch.randn_like(t, device=self.device) * 0.1,
            },
        )

        @torch.compile(dynamic=True)
        def test_model(mat1, mat2):
            return guard_safety_op(mat1, mat2)

        torch._dynamo.reset()
        counters.clear()

        # First call: m=8 (divisible by 4) triggers compilation
        mat1 = torch.randn(8, 64, device=self.device, dtype=self.dtype)
        mat2 = torch.randn(64, 32, device=self.device, dtype=self.dtype)
        with config.patch(max_autotune=True, fx_graph_cache=False):
            result = test_model(mat1, mat2)
        torch.testing.assert_close(result, mat1 @ mat2, rtol=1e-1, atol=1e-1)

        # Verify the guard safety counter was incremented
        self.assertGreater(
            counters["inductor"]["custom_op_decomp_guard_skips"],
            0,
            "Expected custom_op_decomp_guard_skips counter to be incremented",
        )

        # Second call: m=7 (NOT divisible by 4). If the unsafe impl was kept
        # and its Mod(m,4)==0 guard leaked, this would crash on the reshape
        # (can't reshape [7, 64] into [4, m_parts, 64] when 7 isn't divisible by 4).
        mat1_odd = torch.randn(7, 64, device=self.device, dtype=self.dtype)
        mat2_odd = torch.randn(64, 32, device=self.device, dtype=self.dtype)
        with config.patch(max_autotune=True, fx_graph_cache=False):
            result_odd = test_model(mat1_odd, mat2_odd)
        torch.testing.assert_close(
            result_odd, mat1_odd @ mat2_odd, rtol=1e-1, atol=1e-1
        )