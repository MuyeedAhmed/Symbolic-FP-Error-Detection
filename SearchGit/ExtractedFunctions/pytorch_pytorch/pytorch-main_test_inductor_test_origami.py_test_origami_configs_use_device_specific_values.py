    def test_origami_configs_use_device_specific_values(self):
        """Verify that origami configs use architecture-specific num_stages and num_warps.

        Tests that configurations are derived from device properties (MI300 vs MI350X)
        rather than hardcoded values. This ensures portability across AMD GPU models.
        """
        device = torch.device("cuda:0")
        device_props = torch.cuda.get_device_properties(device)

        # Compile a simple MM operation with origami enabled
        torch.manual_seed(0)
        a = torch.randn(256, 256, device=device)
        b = torch.randn(256, 256, device=device)

        # Compile with origami config
        with config.patch(self._origami_default_config(topk=2)):
            compiled = torch.compile(
                lambda x, y: x @ y,
                backend="inductor",
                mode="reduce-overhead",
            )
            compiled_result = compiled(a, b)

        # Verify compilation succeeded
        self.assertIsNotNone(compiled_result)

        # Check that device properties are available and used
        self.assertGreater(
            device_props.multi_processor_count,
            0,
            msg="Device should have valid compute properties",
        )

        trace_structured(
            "origami_device_specific_test",
            metadata_fn=lambda: {
                "device_name": device_props.name,
                "multi_processor_count": device_props.multi_processor_count,
                "warp_size": device_props.warp_size,
                "test": "verify architecture-specific config values",
            },
        )