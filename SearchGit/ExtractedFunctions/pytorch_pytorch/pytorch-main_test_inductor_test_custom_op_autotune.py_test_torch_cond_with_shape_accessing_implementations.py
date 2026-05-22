    def test_torch_cond_with_shape_accessing_implementations(self):
        """Test torch.cond dispatch with implementations that access tensor shapes.

        Validates that implementations like decompose_k that access tensor shapes
        (e.g., `m, k = mat1.shape`) work correctly with torch.cond dispatch.
        The fix uses _build_cond_dispatch_graph to pre-trace each implementation.
        """
        test_op_name = f"test_lib::shape_access_cond_{id(self)}"

        def shape_accessing_impl(
            mat1: torch.Tensor, mat2: torch.Tensor
        ) -> torch.Tensor:
            m, k = mat1.shape  # Shape access that would break naive make_fx
            n = mat2.shape[1]
            k_splits = 4
            if k % k_splits == 0:
                k_parts = k // k_splits
                a = torch.permute(mat1.reshape(m, k_splits, k_parts), (1, 0, 2))
                b = mat2.reshape(k_splits, k_parts, n)
                return torch.sum(torch.bmm(a, b), dim=0)
            return mat1 @ mat2

        def simple_impl(mat1: torch.Tensor, mat2: torch.Tensor) -> torch.Tensor:
            return mat1 @ mat2

        @torch.library.custom_op(test_op_name, mutates_args=())
        def shape_access_op(mat1: torch.Tensor, mat2: torch.Tensor) -> torch.Tensor:
            return mat1 @ mat2

        @shape_access_op.register_fake
        def _(mat1: torch.Tensor, mat2: torch.Tensor):
            return torch.empty(
                mat1.shape[0], mat2.shape[1], device=mat1.device, dtype=mat1.dtype
            )

        register_custom_op_autotuning(
            shape_access_op,
            configs=[CustomOpConfig(simple_impl), CustomOpConfig(shape_accessing_impl)],
            name="shape_access_autotuned",
            dispatch_on={"tensor_name": "mat1", "dim": 0},
            split_points=[4, 16],
            input_gen_fns={
                "mat1": lambda t: torch.randn_like(t, device=self.device) * 0.1,
                "mat2": lambda t: torch.randn_like(t, device=self.device) * 0.1,
            },
        )

        test_mat1 = torch.randn(8, 64, device=self.device, dtype=self.dtype)
        test_mat2 = torch.randn(64, 32, device=self.device, dtype=self.dtype)

        @torch.compile(dynamic=True)
        def test_model(mat1, mat2):
            return shape_access_op(mat1, mat2)

        torch._dynamo.reset()
        with config.patch(max_autotune=True, fx_graph_cache=False):
            result = test_model(test_mat1, test_mat2)

        torch.testing.assert_close(result, test_mat1 @ test_mat2, rtol=1e-1, atol=1e-1)