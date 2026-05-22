    def _solve_eigen_covariance(
        self, alpha, y, sqrt_sw, eigvals, V, X, X_mean, XT_y, XT_sqrt_sw
    ):
        """Compute looe and coef when we have a decomposition of X' X"""
        D = 1 / (eigvals + alpha)
        Hinv = (V * D) @ V.T
        Hinv_XT_y = Hinv @ XT_y
        Hinv_XT_sqrt_sw = Hinv @ XT_sqrt_sw
        X_Hinv_XT_y = safe_sparse_dot(X, Hinv_XT_y, dense_output=True)
        X_Hinv_XT_sqrt_sw = safe_sparse_dot(X, Hinv_XT_sqrt_sw, dense_output=True)
        if self.fit_intercept and sparse.issparse(X):
            # centered = X - sqrt_sw X_mean'
            if y.ndim == 2:
                X_Hinv_XT_y -= sqrt_sw[:, None] * (X_mean @ Hinv_XT_y)
            else:
                X_Hinv_XT_y -= sqrt_sw * (X_mean @ Hinv_XT_y)
            X_Hinv_XT_sqrt_sw -= sqrt_sw * (X_mean @ Hinv_XT_sqrt_sw)
        alpha_c = y - X_Hinv_XT_y
        alpha_d = 1 - self._sparse_multidot_diag(X, Hinv, X_mean, sqrt_sw)
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            alpha_Ginv_sqrt_sw = sqrt_sw - X_Hinv_XT_sqrt_sw
            alpha_d -= alpha_Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            alpha_d = alpha_d[:, None]
        looe = alpha_c / alpha_d
        coef = Hinv_XT_y
        return looe, coef