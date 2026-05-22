    def test_numerical_check_python_binding_tunableop(self, device, dtype):
        with self._tunableop_ctx():
            torch.cuda.tunable.enable(True)
            torch.cuda.tunable.set_numerical_check_tolerances(True)

            a = torch.randn(128, 128, device='cuda')
            b = torch.randn(128, 128, device='cuda')

            _ = a @ b

        with self._tunableop_ctx():
            torch.cuda.tunable.enable(True)
            with self.assertRaisesRegex(RuntimeError, r"positive"):
                torch.cuda.tunable.set_numerical_check_tolerances(True, -1e-5, 1e5)
            with self.assertRaisesRegex(RuntimeError, r"positive"):
                torch.cuda.tunable.set_numerical_check_tolerances(True, 1e-5, -1e5)
            with self.assertRaisesRegex(RuntimeError, r"positive"):
                torch.cuda.tunable.set_numerical_check_tolerances(True, -1e-5, -1e5)