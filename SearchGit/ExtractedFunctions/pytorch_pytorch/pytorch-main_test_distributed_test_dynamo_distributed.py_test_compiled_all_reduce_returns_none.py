    def test_compiled_all_reduce_returns_none(self):
        def fn(x, w):
            result = dist.all_reduce(x, async_op=False)
            assert result is None  # noqa: S101
            return x @ w

        x = torch.randn(4, 4, device=self.device)
        w = torch.randn(4, 4, device=self.device)
        expected = x @ w
        compiled_fn = torch.compile(fn, backend="aot_eager", fullgraph=True)
        actual = compiled_fn(x, w)
        self.assertEqual(actual, expected)