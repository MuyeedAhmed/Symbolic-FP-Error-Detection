    def test_stack_transform(self):
        x1 = -1 * torch.arange(1, 101, dtype=torch.float)
        x2 = (torch.arange(1, 101, dtype=torch.float) - 1) / 100
        x3 = torch.arange(1, 101, dtype=torch.float)
        t1, t2, t3 = ExpTransform(), AffineTransform(1, 100), identity_transform
        dim = 0
        x = torch.stack([x1, x2, x3], dim=dim)
        t = StackTransform([t1, t2, t3], dim=dim)
        actual_dom_check = t.domain.check(x)
        expected_dom_check = torch.stack(
            [t1.domain.check(x1), t2.domain.check(x2), t3.domain.check(x3)], dim=dim
        )
        self.assertEqual(expected_dom_check, actual_dom_check)
        actual = t(x)
        expected = torch.stack([t1(x1), t2(x2), t3(x3)], dim=dim)
        self.assertEqual(expected, actual)
        y1 = torch.arange(1, 101, dtype=torch.float)
        y2 = torch.arange(1, 101, dtype=torch.float)
        y3 = torch.arange(1, 101, dtype=torch.float)
        y = torch.stack([y1, y2, y3], dim=dim)
        actual_cod_check = t.codomain.check(y)
        expected_cod_check = torch.stack(
            [t1.codomain.check(y1), t2.codomain.check(y2), t3.codomain.check(y3)],
            dim=dim,
        )
        self.assertEqual(actual_cod_check, expected_cod_check)
        actual_inv = t.inv(x)
        expected_inv = torch.stack([t1.inv(x1), t2.inv(x2), t3.inv(x3)], dim=dim)
        self.assertEqual(expected_inv, actual_inv)
        actual_jac = t.log_abs_det_jacobian(x, y)
        expected_jac = torch.stack(
            [
                t1.log_abs_det_jacobian(x1, y1),
                t2.log_abs_det_jacobian(x2, y2),
                t3.log_abs_det_jacobian(x3, y3),
            ],
            dim=dim,
        )
        self.assertEqual(actual_jac, expected_jac)