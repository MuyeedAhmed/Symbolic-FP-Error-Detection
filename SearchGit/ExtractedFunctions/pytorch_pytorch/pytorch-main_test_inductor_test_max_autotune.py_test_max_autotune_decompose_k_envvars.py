    def test_max_autotune_decompose_k_envvars(
        self, num_decompose_k_splits, decompose_k_threshold
    ):
        shapes = [(32, 32, 32768), (32, 32, 256)]
        for M, N, K in shapes:
            get_k_splits.cache_clear()
            use_decompose_k_choice.cache_clear()
            a = torch.randn(M, K, dtype=torch.float16, device=GPU_TYPE)
            b = torch.randn(K, N, dtype=torch.float16, device=GPU_TYPE)

            with config.patch(
                {
                    "triton.num_decompose_k_splits": num_decompose_k_splits,
                    "triton.decompose_k_threshold": decompose_k_threshold,
                }
            ):
                compiled_func = torch.compile(lambda a, b: a @ b)
                _, code = run_and_get_code(compiled_func, a, b)

                decompose_count = 0
                for codegen in code:
                    if "benchmark_decompose_k_mm" in codegen:
                        decompose_count += 1

                if (
                    K // M < decompose_k_threshold
                    or K // N < decompose_k_threshold
                    or num_decompose_k_splits == 0
                ):
                    self.assertEqual(decompose_count, 0)
                else:
                    self.assertTrue(decompose_count > 0)
                    self.assertTrue(decompose_count <= num_decompose_k_splits)