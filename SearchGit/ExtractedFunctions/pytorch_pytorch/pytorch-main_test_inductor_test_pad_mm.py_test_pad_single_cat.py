    def test_pad_single_cat(self):
        @torch.compile()
        def foo(x, y):
            return x @ y

        inps = [torch.rand([5, 5], device=GPU_TYPE) for _ in range(2)]
        out = foo(*inps)
        self.assertEqual(out, inps[0] @ inps[1])