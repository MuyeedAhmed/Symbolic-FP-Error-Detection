    def test_permuted_inputs(self):
        B, M, N = 4, 8, 16
        cases = [
            (
                torch.randn(M, B, 1, device=GPU_TYPE).permute(1, 0, 2),
                torch.randn(B, 1, N, device=GPU_TYPE),
            ),
            (
                torch.randn(B, M, 1, device=GPU_TYPE),
                torch.randn(N, B, 1, device=GPU_TYPE).permute(1, 2, 0),
            ),
            (
                torch.randn(M, B, 1, device=GPU_TYPE).permute(1, 0, 2),
                torch.randn(N, B, 1, device=GPU_TYPE).permute(1, 2, 0),
            ),
        ]
        for a, b in cases:
            self.assertEqual(torch.bmm(a, b), a @ b)