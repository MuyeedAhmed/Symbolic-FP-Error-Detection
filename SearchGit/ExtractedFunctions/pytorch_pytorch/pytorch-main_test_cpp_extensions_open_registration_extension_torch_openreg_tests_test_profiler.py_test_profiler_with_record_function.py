    def test_profiler_with_record_function(self):
        """Test profiler with custom record_function annotations."""
        with autograd_profile(use_device="openreg") as prof:
            with record_function("openreg_custom_operation"):
                x = torch.randn(10, 10, device="openreg")
                y = torch.randn(10, 10, device="openreg")
                x @ y

        events = prof.function_events
        event_names = [e.name for e in events]
        self.assertTrue(any("openreg_custom_operation" in name for name in event_names))