    def test_cat_transform(self):
        x1 = -1 * torch.arange(1, 101, dtype=torch.float).view(-1, 100)
        x2 = (torch.arange(1, 101, dtype=torch.float).view(-1, 100) - 1) / 100
        x3 = torch.arange(1, 101, dtype=torch.float).view(-1, 100)
        t1, t2, t3 = ExpTransform(), AffineTransform(1, 100), identity_transform
        dim = 0
        x = torch.cat([x1, x2, x3], dim=dim)
        t = CatTransform([t1, t2, t3], dim=dim)
        actual_dom_check = t.domain.check(x)
        expected_dom_check = torch.cat(
            [t1.domain.check(x1), t2.domain.check(x2), t3.domain.check(x3)], dim=dim
        )
        self.assertEqual(expected_dom_check, actual_dom_check)
        actual = t(x)
        expected = torch.cat([t1(x1), t2(x2), t3(x3)], dim=dim)
        self.assertEqual(expected, actual)
        y1 = torch.arange(1, 101, dtype=torch.float).view(-1, 100)
        y2 = torch.arange(1, 101, dtype=torch.float).view(-1, 100)
        y3 = torch.arange(1, 101, dtype=torch.float).view(-1, 100)
        y = torch.cat([y1, y2, y3], dim=dim)
        actual_cod_check = t.codomain.check(y)
        expected_cod_check = torch.cat(
            [t1.codomain.check(y1), t2.codomain.check(y2), t3.codomain.check(y3)],
            dim=dim,
        )
        self.assertEqual(actual_cod_check, expected_cod_check)
        actual_inv = t.inv(y)
        expected_inv = torch.cat([t1.inv(y1), t2.inv(y2), t3.inv(y3)], dim=dim)
        self.assertEqual(expected_inv, actual_inv)
        actual_jac = t.log_abs_det_jacobian(x, y)
        expected_jac = torch.cat(
            [
                t1.log_abs_det_jacobian(x1, y1),
                t2.log_abs_det_jacobian(x2, y2),
                t3.log_abs_det_jacobian(x3, y3),
            ],
            dim=dim,
        )
        self.assertEqual(actual_jac, expected_jac)