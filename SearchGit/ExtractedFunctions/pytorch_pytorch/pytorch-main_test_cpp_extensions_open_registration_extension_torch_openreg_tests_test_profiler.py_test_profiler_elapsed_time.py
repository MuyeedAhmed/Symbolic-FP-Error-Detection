    def test_profiler_elapsed_time(self):
        """Test that event elapsed time is correctly calculated."""
        with autograd_profile(use_device="openreg") as prof:
            stream = torch.Stream(device="openreg")

            with stream:
                x = torch.randn(100, 100, device="openreg")
                y = torch.randn(100, 100, device="openreg")
                # Multiple matmuls to ensure measurable time
                for _ in range(10):
                    z = x @ y
                    x = z

            stream.synchronize()

        events = prof.function_events
        # Check that operations have non-zero duration
        compute_events = [
            e for e in events if "aten::mm" in e.name or "aten::matmul" in e.name
        ]
        if compute_events:
            total_time = sum(e.cpu_time_total for e in compute_events)
            self.assertGreater(total_time, 0)