    def test_no_recursion_in_aten_backward_formula(self):
        """Backward of bmm calls bmm again; the fallback must bypass the router.

        If the router weren't bypassed, aten's bmm backward formula would
        re-enter it on every autograd step, either recursing forever or
        silently calling the override impl in the grad path.
        """

        def cond(*a, **k):
            # False so every call hits the fallback — this isolates the
            # "native kernel called via fallback must not re-enter the router"
            # property.
            return False

        def impl(a, b):
            raise AssertionError("impl should not be called when cond=False")

        self.registry.register_op_override("test_dsl", "aten", "bmm", "CPU", cond, impl)
        self._install("bmm", "CPU")

        a = torch.randn(2, 3, 4, requires_grad=True)
        b = torch.randn(2, 4, 5, requires_grad=True)
        expected = a.detach() @ b.detach()
        out = torch.bmm(a, b)
        self.assertEqual(out, expected)

        out.sum().backward()
        self.assertEqual(a.grad.shape, a.shape)
        self.assertEqual(b.grad.shape, b.shape)
        self.assertTrue((a.grad != 0).any())
        self.assertTrue((b.grad != 0).any())