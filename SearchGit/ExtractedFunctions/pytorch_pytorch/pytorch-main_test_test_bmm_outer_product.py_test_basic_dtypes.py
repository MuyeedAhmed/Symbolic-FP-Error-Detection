    def test_basic_dtypes(self):
        for dtype in [torch.float32, torch.float16, torch.bfloat16]:
            with self.subTest(dtype=dtype):
                a = torch.randn(4, 8, 1, device=GPU_TYPE, dtype=dtype)
                b = torch.randn(4, 1, 16, device=GPU_TYPE, dtype=dtype)
                self.assertEqual(torch.bmm(a, b), a @ b)