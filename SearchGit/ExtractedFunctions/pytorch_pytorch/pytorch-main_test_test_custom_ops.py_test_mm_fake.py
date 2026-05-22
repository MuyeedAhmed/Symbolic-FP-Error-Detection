    def test_mm_fake(self):
        with torch._subclasses.fake_tensor.FakeTensorMode():
            x = torch.randn(2, 3, requires_grad=True, device="cpu")
            y = torch.randn(3, 5, device="cpu")
            result = torch.ops.aten.mm.default(x, y)
            self.assertEqual(result.shape, (x @ y).shape)