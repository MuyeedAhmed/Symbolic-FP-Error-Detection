    def test_compiled_all_gather_returns_none(self):
        def fn(tensor_list, tensor, w):
            result = dist.all_gather(tensor_list, tensor, async_op=False)
            assert result is None  # noqa: S101
            return tensor_list[0] @ w

        tensor = torch.randn(4, 4, device=self.device)
        tensor_list = [torch.empty(4, 4, device=self.device)]
        w = torch.randn(4, 4, device=self.device)
        compiled_fn = torch.compile(fn, backend="aot_eager", fullgraph=True)
        actual = compiled_fn(tensor_list, tensor, w)
        self.assertEqual(actual, tensor @ w)