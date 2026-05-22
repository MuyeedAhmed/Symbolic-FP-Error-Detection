    def test_prologue(self):
        if not is_big_gpu():
            raise unittest.SkipTest("requires large gpu to max-autotune")

        def gn(x, y):
            return torch.mul(x, y) + (y - 1)

        def fn(x, y, z):
            return (
                invoke_quant_tracer(
                    gn, x, y, scheme="nf4", quant_options=invoke_quant_tracer
                )
                @ z
            )

        x = torch.randn(
            64, 64, requires_grad=False, device=GPU_TYPE, dtype=torch.float16
        )
        # make this a no-op to ensure equivalent numerics
        y = torch.randn(
            64, 64, requires_grad=False, device=GPU_TYPE, dtype=torch.float16
        ).fill_(1.0)
        z = torch.randn(
            64, 64, requires_grad=False, device=GPU_TYPE, dtype=torch.float16
        )
        ref = gn(x, y) @ z

        x_clone = x.clone().detach().requires_grad_(False)
        y_clone = y.clone().detach().requires_grad_(False)
        z_clone = z.clone().detach().requires_grad_(False)
        torch._dynamo.reset()
        with torch.no_grad(), config.patch(max_autotune_gemm_backends="TRITON"):
            fn_c = torch.compile(fn, mode="max-autotune-no-cudagraphs")
            res, code = run_and_get_code(fn_c, x_clone, y_clone, z_clone)

            FileCheck().check("k_idx in range").check_not("tl.float32").check(
                "tl.dot"
            ).run(code[0])

            self.assertEqual(ref, res)