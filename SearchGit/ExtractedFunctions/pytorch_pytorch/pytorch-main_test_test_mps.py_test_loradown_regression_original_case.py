    def test_loradown_regression_original_case(self):
        a = torch.rand(2, 1025, device='mps', dtype=torch.half)
        b = torch.rand(2, 1041, device='mps', dtype=torch.half)[:, :1025].t()
        result = a @ b
        self.assertEqual(result.shape, (2, 2))

        self.assertFalse(torch.isnan(result).any())
        self.assertFalse(torch.isinf(result).any())