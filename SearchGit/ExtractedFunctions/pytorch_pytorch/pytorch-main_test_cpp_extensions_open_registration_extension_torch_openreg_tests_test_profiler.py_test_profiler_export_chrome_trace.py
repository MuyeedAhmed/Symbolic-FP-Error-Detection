    def test_profiler_export_chrome_trace(self):
        """Test exporting profiler results to Chrome trace format."""
        with autograd_profile(use_device="openreg") as prof:
            x = torch.randn(10, 10, device="openreg")
            y = x + x
            x @ y

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            trace_file = f.name

        try:
            prof.export_chrome_trace(trace_file)

            # Verify the trace file is valid JSON
            with open(trace_file) as f:
                trace_data = json.load(f)

            if isinstance(trace_data, dict):
                self.assertIn("traceEvents", trace_data)
                events = trace_data["traceEvents"]
            else:
                events = trace_data
            self.assertGreater(len(events), 0)
        finally:
            import os

            if os.path.exists(trace_file):
                os.remove(trace_file)