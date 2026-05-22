    def test_mm_meta(self):
        x = torch.randn(2, 3, requires_grad=True, device="meta")
        y = torch.randn(3, 5, device="meta")
        result = torch.ops.aten.mm.default(x, y)
        self.assertEqual(result.shape, (x @ y).shape)