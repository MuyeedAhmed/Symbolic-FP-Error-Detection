    def test_profiler_matmul_dim_fp16_pattern(self):
        cases = (
            (1, torch.randn((201, 201), device="cuda", dtype=torch.float16)),
            (1, torch.randn((3, 97, 97), device="cuda", dtype=torch.float16)),
            (0, torch.randn((200, 200), device="cuda", dtype=torch.float16)),
            (0, torch.randn((3, 200, 200), device="cuda", dtype=torch.float16)),
        )
        num_matched = []
        for _, x in cases:
            with profile(with_stack=True, record_shapes=True) as prof:
                x @ x
            pattern = MatMulDimInFP16Pattern(prof)
            num_matched.append(len(pattern.matched_events()))
        self.assertEqual(num_matched, [i for i, _ in cases])