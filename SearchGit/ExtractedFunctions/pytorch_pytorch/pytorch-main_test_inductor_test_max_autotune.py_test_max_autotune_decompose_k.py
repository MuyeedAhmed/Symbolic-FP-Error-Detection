    def test_max_autotune_decompose_k(self, sizes, dtype, dynamic):
        # UT specific change to force testing decompose K feature on ROCm until
        # enabled by default, same strategy as #169948
        with config.patch(_DECOMPOSE_K_PATCH_ROCM):
            fp16_red_setting = (
                torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction
            )
            bf16_red_setting = (
                torch.backends.cuda.matmul.allow_bf16_reduced_precision_reduction
            )
            torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction = False
            torch.backends.cuda.matmul.allow_bf16_reduced_precision_reduction = False

            M, N, K = sizes

            atol = 1e-4
            rtol = 1e-4
            # K can be huge huge, this is why the data distribution is set to iid N(0, K ** 0.5),
            # which makes the result of reductions distributed as N(0, 1).
            a, b = self._make_matrices(
                M,
                K,
                N,
                dtype=dtype,
                device=GPU_TYPE,
                requires_grad=True,
            )

            possible_splits = range(2, min(K // M, K // N) + 1)

            divisors = {split for split in possible_splits if K % split == 0}

            def check_divisors(code):
                for kernel in code:
                    if "decompose_k" in kernel:
                        divisor_found = False
                        for divisor in divisors:
                            if f"{divisor}_split" in kernel:
                                divisor_found = True
                                break

                        self.assertTrue(
                            divisor_found,
                            f"Could not find a split in {divisors} in {kernel}",
                        )

            compiled_func = torch.compile(lambda a, b: a @ b, dynamic=dynamic)
            # We assume with the large k dim relative to m, n, decompose_k will be most performant
            out, code = run_and_get_code(compiled_func, a, b)

            if dynamic:
                FileCheck().check_not("extern_kernels.bmm_dtype").check_not(
                    "decompose_k"
                ).run(code[0])
            else:
                FileCheck().check("extern_kernels.bmm_dtype").check_regex(
                    "triton_.*_fused_.*.run"
                ).check("decompose_k").run(code[0])
                check_divisors(code)
                torch.testing.assert_close(out, a @ b, atol=atol, rtol=rtol)

            # Test adding epilogue also equivalent to eager
            compiled_func = torch.compile(lambda a, b: (a @ b).relu(), dynamic=dynamic)
            out, code = run_and_get_code(compiled_func, a, b)
            if dynamic:
                FileCheck().check_not("extern_kernels.bmm_dtype").check_not(
                    "decompose_k"
                ).run(code[0])
            else:
                FileCheck().check("extern_kernels.bmm_dtype").check_regex(
                    "triton_.*_fused_.*.run"
                ).check("decompose_k").run(code[0])
                check_divisors(code)
                torch.testing.assert_close(
                    compiled_func(a, b), (a @ b).relu(), atol=atol, rtol=rtol
                )

            # Test adding reinterpret view before subgraph
            a = a.transpose(0, 1)
            compiled_func = torch.compile(
                lambda a, b: (a.transpose(0, 1) @ b).relu(), dynamic=dynamic
            )
            out, code = run_and_get_code(compiled_func, a, b)

            if dynamic:
                FileCheck().check_not("extern_kernels.bmm_dtype").check_not(
                    "decompose_k"
                ).run(code[0])
            else:
                FileCheck().check("extern_kernels.bmm_dtype").check_regex(
                    "triton_.*_fused_.*_0.run"
                ).check("decompose_k").run(code[0])
                check_divisors(code)
                torch.testing.assert_close(
                    compiled_func(a, b),
                    (a.transpose(0, 1) @ b).relu(),
                    atol=atol,
                    rtol=rtol,
                )

            torch.backends.cuda.matmul.allow_fp16_reduced_precision_reduction = (
                fp16_red_setting
            )
            torch.backends.cuda.matmul.allow_bf16_reduced_precision_reduction = (
                bf16_red_setting
            )