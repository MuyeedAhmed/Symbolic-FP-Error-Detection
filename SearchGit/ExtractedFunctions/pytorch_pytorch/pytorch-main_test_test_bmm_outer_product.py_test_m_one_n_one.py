    def test_m_one_n_one(self):
        a = torch.randn(8, 1, 1, device=GPU_TYPE)
        b = torch.randn(8, 1, 1, device=GPU_TYPE)
        self.assertEqual(torch.bmm(a, b), a @ b)