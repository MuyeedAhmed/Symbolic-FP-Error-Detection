    def test_profiler_fp32_matmul_pattern(self):
        x = torch.ones((100, 100), device="cuda")
        with profile(with_stack=True) as prof:
            x = x @ x
        pattern = FP32MatMulPattern(prof)
        has_tf32 = 0 if pattern.skip else 1
        num_matched = len(pattern.matched_events())
        self.assertEqual(num_matched, has_tf32)