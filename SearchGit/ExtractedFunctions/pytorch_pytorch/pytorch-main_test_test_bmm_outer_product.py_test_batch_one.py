    def test_batch_one(self):
        a = torch.randn(1, 64, 1, device=GPU_TYPE)
        b = torch.randn(1, 1, 128, device=GPU_TYPE)
        self.assertEqual(torch.bmm(a, b), a @ b)