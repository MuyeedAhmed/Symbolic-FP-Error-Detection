    def _solve_eigen_gram(
        self, alpha, y, sqrt_sw, eigvals, Q, QT_y, QT_sqrt_sw, XT, X_mean
    ):
        """Compute looe and coef when we have a decomposition of X X'"""
        D = 1.0 / (eigvals + alpha)
        c = Q @ self._diag_dot(D, QT_y)
        d = self._decomp_diag(D, Q)
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            Ginv_sqrt_sw = Q @ self._diag_dot(D, QT_sqrt_sw)
            d -= Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            d = d[:, None]
        XT_c = XT @ c
        if self.fit_intercept and sparse.issparse(XT):
            # centered matrix = X - sqrt_sw X_mean'
            if y.ndim == 2:
                XT_c -= X_mean[:, None] * (sqrt_sw @ c)
            else:
                XT_c -= X_mean * (sqrt_sw @ c)
        looe = c / d
        coef = XT_c
        return looe, coef