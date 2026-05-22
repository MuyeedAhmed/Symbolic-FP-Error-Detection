    def test_profiler_with_data_transfer(self):
        """Test profiler with CPU-OpenReg data transfers."""
        with autograd_profile(use_device="openreg") as prof:
            # CPU tensor
            x_cpu = torch.randn(50, 50)

            # Transfer to OpenReg
            with record_function("to_openreg"):
                x_openreg = x_cpu.to("openreg")

            # Compute on OpenReg
            with record_function("compute_openreg"):
                y = x_openreg @ x_openreg

            # Transfer back to CPU
            with record_function("to_cpu"):
                y.cpu()

        events = prof.function_events
        event_names = [e.name for e in events]

        self.assertTrue(any("to_openreg" in name for name in event_names))
        self.assertTrue(any("compute_openreg" in name for name in event_names))
        self.assertTrue(any("to_cpu" in name for name in event_names))