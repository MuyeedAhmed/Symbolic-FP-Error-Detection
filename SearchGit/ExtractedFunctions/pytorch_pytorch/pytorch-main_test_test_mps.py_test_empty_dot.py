    def test_empty_dot(self):
        # just to check that it doesn't crash
        a = torch.rand((0), device="mps")
        b = torch.rand((0), device="mps")
        self.assertEqual(a.dot(b), a.cpu().dot(b.cpu()))