        def test_amp_cache_disabled(self):
            @torch.compile()
            def foo(x):
                return x + x

            for _ in range(3):
                out = foo(torch.rand([4, 4], device="cuda", requires_grad=True))

            # amp cache for cudagraph outputs should be disabled
            t2 = torch.rand([4, 4], device="cuda")

            with torch.cuda.amp.autocast():
                run_once = out @ t2

                out.detach().zero_()

                run_twice = out @ t2

                self.assertNotEqual(run_once, run_twice)