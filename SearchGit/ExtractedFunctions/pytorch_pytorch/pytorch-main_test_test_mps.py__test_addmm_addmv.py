    def _test_addmm_addmv(self, f, t, m, v, *, alpha=None, beta=None, transpose_out=False):
        dtype = t.dtype
        numpy_dtype = dtype
        alpha = 1.2 if alpha is None else alpha
        beta = 0.8 if beta is None else beta
        res1 = f(t, m, v, alpha=alpha, beta=beta)
        res2 = torch.full_like(res1, math.nan)
        if transpose_out:
            res2 = res2.t().clone(memory_format=torch.contiguous_format).t()
        f(t, m, v, alpha=alpha, beta=beta, out=res2)
        res3 = alpha * (m.to(numpy_dtype).cpu().numpy() @ v.to(numpy_dtype).cpu().numpy())
        if beta != 0:
            res3 += (torch.mul(t, beta)).to(numpy_dtype).cpu().numpy()
        res3 = torch.from_numpy(res3).to(dtype)
        self.assertEqual(res1, res2)
        self.assertEqual(res1, res3)