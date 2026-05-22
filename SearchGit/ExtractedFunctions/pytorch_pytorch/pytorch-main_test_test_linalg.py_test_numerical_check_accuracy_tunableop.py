    def test_numerical_check_accuracy_tunableop(self, device, dtype):
        shapes = [(127, 193, 61), (251, 317, 73), (89, 149, 41)]
        atol, rtol = 1e-2, 1e-1

        for (m, k, n) in shapes:
            a = torch.randn(m, k, device='cuda')
            b = torch.randn(k, n, device='cuda')
            torch.cuda.tunable.enable(False)
            torch.cuda.tunable.set_numerical_check_tolerances(False)
            C_baseline = a @ b
            with self._tunableop_ctx():
                torch.cuda.tunable.enable(True)
                torch.cuda.tunable.set_numerical_check_tolerances(True, atol, rtol)
                C_numeric = a @ b
            self.assertTrue(torch.allclose(C_baseline, C_numeric, atol=atol, rtol=rtol))