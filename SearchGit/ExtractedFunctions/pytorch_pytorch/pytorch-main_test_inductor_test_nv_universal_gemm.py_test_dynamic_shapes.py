    def test_dynamic_shapes(self):
        """Stress test dynamic shapes with extreme variations."""

        def matmul(a, b):
            return a @ b

        torch._dynamo.reset()

        with config.patch(_nvgemm_config(nvgemm_max_profiling_configs=2)):
            compiled_fn = torch.compile(matmul, dynamic=True)

            shapes = [
                (4, 4, 4, False),
                (16, 16, 16, True),
                (2048, 64, 128, True),
                (4, 4, 4, False),  # Unsupported again
                (64, 2048, 128, True),
                (128, 128, 2048, True),
                (2048, 2048, 512, True),
                (16, 16, 16, True),
            ]

            for m, n, k, supported in shapes:
                a = torch.randn(m, k, dtype=torch.bfloat16, device="cuda")
                b = torch.randn(k, n, dtype=torch.bfloat16, device="cuda")
                if not supported:
                    with self.assertRaisesRegex(
                        Exception, "NoValidChoicesError|no valid choice"
                    ):
                        compiled_fn(a, b)
                else:
                    result = compiled_fn(a, b)
                    torch.testing.assert_close(result, a @ b)