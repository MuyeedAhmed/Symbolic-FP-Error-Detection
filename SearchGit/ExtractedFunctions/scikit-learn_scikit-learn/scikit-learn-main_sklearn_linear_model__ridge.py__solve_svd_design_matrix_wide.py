    def _solve_svd_design_matrix_wide(
        self, alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
    ):
        """Compute looe and coef when we have an SVD decomposition of X.

        Wide X case (n_samples < n_features).
        """
        alpha_D = alpha / (singvals**2 + alpha)
        alpha_c = U @ self._diag_dot(alpha_D, UT_y)
        alpha_d = self._decomp_diag(alpha_D, U)
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            alpha_Ginv_sqrt_sw = U @ self._diag_dot(alpha_D, UT_sqrt_sw)
            alpha_d -= alpha_Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            # handle case where y is 2-d
            alpha_d = alpha_d[:, None]
        looe = alpha_c / alpha_d
        coef = V @ self._diag_dot(singvals / (singvals**2 + alpha), UT_y)
        return looe, coef