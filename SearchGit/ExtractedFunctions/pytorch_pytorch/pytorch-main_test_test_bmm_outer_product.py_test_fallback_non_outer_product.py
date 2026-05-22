    def test_fallback_non_outer_product(self):
        a = torch.randn(4, 8, 16, device=GPU_TYPE)
        b = torch.randn(4, 16, 32, device=GPU_TYPE)
        self.assertEqual(torch.bmm(a, b), a @ b, atol=1e-5, rtol=1.3e-6)