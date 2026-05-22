    def test_mm(self):
        x = torch.randn(2, 3, requires_grad=True)
        y = torch.randn(3, 5)
        result = torch.ops.aten.mm.default(x, y)
        self.assertEqual(result, x @ y)