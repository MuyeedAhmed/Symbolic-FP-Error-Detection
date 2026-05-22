        def test_empty_inputs(lhs_layout, rhs_layout):
            xd = torch.rand(10, 0, device=device, dtype=dtype)
            yd = xd.transpose(-2, -1)
            zd = torch.rand(0, 0, device=device, dtype=dtype)

            xls, yls, zls = (t.to_sparse(layout=lhs_layout) for t in (xd, yd, zd))
            xrs, yrs, zrs = (t.to_sparse(layout=rhs_layout) for t in (xd, yd, zd))

            for ls, rs, ld, rd in [(xls, yrs, xd, yd), (xls, zrs, xd, zd), (zls, yrs, zd, yd), (zls, zrs, zd, zd)]:
                res_sparse = ls @ rs
                res_dense = ld @ rd
                self.assertEqual(res_sparse.to_dense(), res_dense)