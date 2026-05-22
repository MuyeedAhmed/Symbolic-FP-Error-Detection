    def test_default_dispatching(self):
        # check that the correct implementation is dispatched,
        # and ops do not modify inputs when using the default overload
        w = torch.randn(3, 3)
        x = torch.randn(2, 3)
        x1 = x.clone()

        with _custom_mm2.set_priority(["regular"]):
            result_regular = _custom_mm2(x, w)

        with _custom_mm2.set_priority(["inplace"]):
            result_inplace = _custom_mm2(x, w)

        # check that x was not modified by either impl
        torch.testing.assert_close(x, x1, atol=0, rtol=0)

        torch.testing.assert_close(result_inplace, x1 @ w + 2)
        torch.testing.assert_close(result_regular, x1 @ w + 1)