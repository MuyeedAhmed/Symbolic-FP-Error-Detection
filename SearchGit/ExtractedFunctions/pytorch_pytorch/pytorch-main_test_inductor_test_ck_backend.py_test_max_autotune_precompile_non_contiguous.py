    def test_max_autotune_precompile_non_contiguous(self, max_autotune_gemm_backends):
        """
        Make sure the matmul with non-contiguous inputs can fallback
        """

        tensor_options = {"device": "cuda", "dtype": torch.float16}

        a = torch.empty_strided((50257, 32768), (1, 50304), **tensor_options)
        b = torch.empty_strided((32768, 768), (768, 1), **tensor_options)

        if "rocm" not in dir(config):
            raise AssertionError("'rocm' not found in dir(config)")

        with (
            config.patch(
                {
                    "max_autotune": True,
                    "autotune_in_subproc": True,
                    "max_autotune_gemm_backends": max_autotune_gemm_backends,
                    "compile_threads": 16,
                    "rocm.ck_dir": self.ck_dir,
                    "rocm.ck_max_profiling_configs": 8,
                    "rocm.ck_tile_max_profiling_configs": 8,
                }
            ),
            tf32_off(),
        ):

            @torch.compile(dynamic=False)
            def mm(a, b):
                return a @ b

            Y_compiled = mm(a, b)
            Y_eager = a @ b
            torch.testing.assert_close(Y_compiled, Y_eager, equal_nan=True)