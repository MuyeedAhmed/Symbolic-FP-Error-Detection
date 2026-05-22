    def test_non_contiguous_input_mm_plus_mm(self):
        x1 = rand_strided((50257, 2048), (1, 50304), device=GPU_TYPE)
        y1 = rand_strided((2048, 768), (768, 1), device=GPU_TYPE)

        x2 = rand_strided((50257, 2048), (1, 50304), device=GPU_TYPE)
        y2 = rand_strided((2048, 768), (768, 1), device=GPU_TYPE)

        @torch.compile(mode="max-autotune")
        def f(x1, y1, x2, y2):
            return x1 @ y1 + x2 @ y2

        ref = x1 @ y1 + x2 @ y2
        act = f(x1, y1, x2, y2)
        torch.testing.assert_close(act, ref, atol=1e-1, rtol=1e-2)