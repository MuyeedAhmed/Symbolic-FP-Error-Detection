    def test_profiler_with_backward_pass(self):
        """Test profiler with autograd operations."""
        with autograd_profile(use_device="openreg", record_shapes=True) as prof:
            x = torch.randn(10, 10, device="openreg", requires_grad=True)
            y = torch.randn(10, 10, device="openreg", requires_grad=True)
            z = (x @ y).sum()
            z.backward()

        events = prof.function_events
        event_names = [e.name for e in events]

        # Check for forward and backward operations
        self.assertTrue(
            any("aten::mm" in name or "aten::matmul" in name for name in event_names)
        )
        self.assertTrue(any("backward" in name.lower() for name in event_names))