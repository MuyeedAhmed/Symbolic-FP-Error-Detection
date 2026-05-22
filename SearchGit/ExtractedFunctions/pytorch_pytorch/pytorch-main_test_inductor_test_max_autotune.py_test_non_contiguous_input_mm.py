    def test_non_contiguous_input_mm(self):
        """
        Make sure the triton template can work with non-contiguous inputs without crash.
        Check https://github.com/pytorch/pytorch/issues/125437 for more details.
        """
        x = rand_strided(
            (50257, 2048), (1, 50304), dtype=torch.bfloat16, device=GPU_TYPE
        )
        y = rand_strided((2048, 768), (768, 1), dtype=torch.bfloat16, device=GPU_TYPE)

        @torch.compile(mode="max-autotune")
        def f(x, y):
            return x @ y

        ref = x @ y
        act = f(x, y)
        torch.testing.assert_close(act, ref, atol=2e-2, rtol=1e-2)