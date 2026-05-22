    def test_compiled_reduce_scatter_tensor_returns_none(self):
        def fn(output, input, w):
            result = dist.reduce_scatter_tensor(output, input, async_op=False)
            assert result is None  # noqa: S101
            return output @ w

        input = torch.randn(4, 4, device=self.device)
        output = torch.empty(4, 4, device=self.device)
        w = torch.randn(4, 4, device=self.device)
        compiled_fn = torch.compile(fn, backend="aot_eager", fullgraph=True)
        actual = compiled_fn(output, input, w)
        self.assertEqual(actual, input @ w)