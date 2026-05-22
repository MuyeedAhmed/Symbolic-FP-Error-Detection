    def test_profiler_with_shapes(self):
        """Test profiler with shape recording enabled."""
        with autograd_profile(use_device="openreg", record_shapes=True) as prof:
            x = torch.randn(10, 20, device="openreg")
            y = torch.randn(20, 30, device="openreg")
            x @ y

        key_averages = prof.key_averages(group_by_input_shape=True)
        self.assertGreater(len(key_averages), 0)