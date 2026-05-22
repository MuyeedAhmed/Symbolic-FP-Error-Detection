    def test_cpu_outer_product_fallback(self):
        a = torch.randn(4, 8, 1)
        b = torch.randn(4, 1, 16)
        self.assertTrue(a.device.type == "cpu")
        self.assertEqual(torch.bmm(a, b), a @ b)