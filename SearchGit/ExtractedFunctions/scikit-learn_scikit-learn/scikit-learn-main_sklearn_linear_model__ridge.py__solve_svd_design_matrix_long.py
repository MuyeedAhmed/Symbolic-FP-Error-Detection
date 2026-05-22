    def _solve_svd_design_matrix_long(
        self, alpha, y, sqrt_sw, singvals, U, V, UT_y, UT_sqrt_sw
    ):
        """Compute looe and coef when we have an SVD decomposition of X.

        Long X case (n_features < n_samples).
        """
        M = alpha / (singvals**2 + alpha) - 1
        alpha_c = U @ self._diag_dot(M, UT_y) + y
        alpha_d = self._decomp_diag(M, U) + 1
        if self.fit_intercept:
            sw_sum = sqrt_sw @ sqrt_sw
            alpha_Ginv_sqrt_sw = U @ self._diag_dot(M, UT_sqrt_sw) + sqrt_sw
            alpha_d -= alpha_Ginv_sqrt_sw * sqrt_sw / sw_sum
        if y.ndim == 2:
            # handle case where y is 2-d
            alpha_d = alpha_d[:, None]
        looe = alpha_c / alpha_d
        coef = V @ self._diag_dot(singvals / (singvals**2 + alpha), UT_y)
        return looe, coef