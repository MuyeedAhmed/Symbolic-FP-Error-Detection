    def test_inplace_dispatching(self):
        # check that the correct implementation is dispatched based on priority,
        # and inplace semantics hold
        w = torch.randn(3, 3)
        x = torch.randn(2, 3)
        x1 = x.clone()

        with _custom_mm2.set_priority(["regular"]):
            result_regular = _custom_mm2.maybe_inplace(x, w)

        # check that the regular op does not modify x
        torch.testing.assert_close(x, x1, atol=0, rtol=0)

        with _custom_mm2.set_priority(["inplace"]):
            result_inplace: Tensor = _custom_mm2.maybe_inplace(x, w)

        # check that the inplace op returns x directly
        assert result_inplace.data_ptr() == x.data_ptr()

        torch.testing.assert_close(result_inplace, x1 @ w + 2)
        torch.testing.assert_close(result_regular, x1 @ w + 1)