    def blas_lapack_ops(self):
        m = torch.randn(3, 3)
        a = torch.randn(10, 3, 4)
        b = torch.randn(10, 4, 3)
        v = torch.randn(3)
        return len(
            torch.addbmm(m, a, b),
            torch.addmm(torch.randn(2, 3), torch.randn(2, 3), torch.randn(3, 3)),
            torch.addmv(torch.randn(2), torch.randn(2, 3), torch.randn(3)),
            torch.addr(torch.zeros(3, 3), v, v),
            torch.baddbmm(m, a, b),
            torch.bmm(a, b),
            torch.chain_matmul(torch.randn(3, 3), torch.randn(3, 3), torch.randn(3, 3)),
            # torch.cholesky(a), # deprecated
            # torch.cholesky_inverse(torch.randn(3, 3)), # had some error
            # torch.cholesky_solve(torch.randn(3, 3), torch.randn(3, 3)),
            torch.dot(v, v),
            # torch.linalg.eig(m), # not build with lapack
            # torch.geqrf(a),
            torch.ger(v, v),
            torch.inner(m, m),
            # torch.inverse(m),
            # torch.det(m),
            # torch.logdet(m),
            # torch.slogdet(m),
            # torch.lstsq(m, m),
            # torch.linalg.lu_factor(m),
            # torch.lu_solve(m, *torch.linalg.lu_factor(m)),
            # torch.lu_unpack(*torch.linalg.lu_factor(m)),
            torch.matmul(m, m),
            torch.matrix_power(m, 2),
            # torch.matrix_rank(m),
            torch.matrix_exp(m),
            torch.mm(m, m),
            torch.mv(m, v),
            # torch.orgqr(a, m),
            # torch.ormqr(a, m, v),
            torch.outer(v, v),
            # torch.pinverse(m),
            # torch.qr(a),
            # torch.solve(m, m),
            # torch.svd(a),
            # torch.svd_lowrank(a),
            # torch.pca_lowrank(a),
            # torch.symeig(a), # deprecated
            # torch.lobpcg(a, b), # not supported
            torch.trapz(m, m),
            torch.trapezoid(m, m),
            torch.cumulative_trapezoid(m, m),
            # torch.triangular_solve(m, m),
            torch.vdot(v, v),
        )