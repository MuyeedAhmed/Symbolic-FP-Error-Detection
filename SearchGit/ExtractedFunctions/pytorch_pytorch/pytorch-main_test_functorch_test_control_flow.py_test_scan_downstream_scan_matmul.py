    def test_scan_downstream_scan_matmul(self, compile_mode, reverse, device, autograd):
        inp = torch.randn(3, 10, 2, device=device, requires_grad=autograd)
        init = torch.randn(3, 2, device=device, requires_grad=autograd)

        for ind in range(2):
            # Chain with matmul
            def chain_fct(inp):
                W = torch.ones(2, 5, device=device)
                o = scan(
                    get_scan_combine_fn("add", False),
                    init,
                    inp,
                    dim=1,
                    reverse=reverse,
                )
                return o[ind] @ W

            fct_cmp = compile_mode_helper(chain_fct, compile_mode)

            expected_result = _fake_scan(
                get_scan_combine_fn("add", False),
                init=init,
                xs=inp,
                dim=1,
                reverse=reverse,
            )[ind] @ torch.ones(2, 5, device=device)
            result = fct_cmp(inp)
            self.assertEqual(result, expected_result)

            if autograd:
                self.check_autograd(result, expected_result, (init, inp))