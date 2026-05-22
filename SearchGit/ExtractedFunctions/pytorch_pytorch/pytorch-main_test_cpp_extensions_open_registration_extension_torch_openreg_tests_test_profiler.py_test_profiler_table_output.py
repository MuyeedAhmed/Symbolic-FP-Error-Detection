    def test_profiler_table_output(self):
        """Test that profiler can generate table output."""
        with autograd_profile(use_device="openreg") as prof:
            x = torch.randn(10, 10, device="openreg")
            y = x + x
            x @ y

        # Should not raise any exceptions
        table = prof.key_averages().table(sort_by="cpu_time_total", row_limit=10)
        self.assertIsInstance(table, str)
        self.assertGreater(len(table), 0)