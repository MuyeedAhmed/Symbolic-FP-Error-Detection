        def test_orthogonal_inputs(lhs_layout, rhs_layout):
            ones = torch.ones(2, 2, device=device, dtype=dtype)
            zeros = torch.zeros(2, 2, device=device, dtype=dtype)
            x = torch.cat((ones, zeros), -1).to_sparse(layout=lhs_layout)
            y = torch.cat((zeros, ones), -2).to_sparse(layout=rhs_layout)
            res = x @ y
            res_expected = torch.zeros(*res.shape, device=device, dtype=dtype, layout=res.layout)
            self.assertEqual(res, res_expected)