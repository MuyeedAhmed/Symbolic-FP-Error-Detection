    def test_compiled_all_to_all_single_returns_none(self):
        def fn(output, input, w):
            result = dist.all_to_all_single(output, input, async_op=False)
            assert result is None  # noqa: S101
            return output @ w

        input = torch.randn(4, 4, device=self.device)
        output = torch.empty(4, 4, device=self.device)
        w = torch.randn(4, 4, device=self.device)
        compiled_fn = torch.compile(fn, backend="aot_eager", fullgraph=True)
        actual = compiled_fn(output, input, w)
        self.assertEqual(actual, input @ w)