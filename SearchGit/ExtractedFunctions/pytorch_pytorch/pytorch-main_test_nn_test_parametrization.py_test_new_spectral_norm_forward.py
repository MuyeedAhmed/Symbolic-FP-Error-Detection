    def test_new_spectral_norm_forward(self):
        input = torch.randn(3, 5)
        m = nn.Linear(5, 7)
        m = torch.nn.utils.parametrizations.spectral_norm(m)
        snm = m.parametrizations.weight[0]
        # naive forward
        _weight = m.parametrizations.weight.original
        _bias, _v = m.bias, snm._v
        _weight_mat = _weight.view(_weight.size(0), -1)
        _u = torch.mv(_weight_mat, _v)
        _u = F.normalize(_u, dim=0, eps=1e-12)
        _v = torch.mv(_weight_mat.t(), _u)
        _v = F.normalize(_v, dim=0, eps=1e-12)
        _weight.data /= torch.dot(_u, torch.matmul(_weight_mat, _v))
        out_hat = torch.nn.functional.linear(input, _weight, _bias)
        expect_out = m(input)
        self.assertEqual(expect_out, out_hat)