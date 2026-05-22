    def _eigen_decompose_covariance(self, X, X_mean, y, sqrt_sw):
        """Eigendecomposition of covariance matrix X' X"""
        xp, is_array_api = get_namespace(X)
        cov = self._compute_covariance(X, X_mean, sqrt_sw)
        eigvals, V = xp.linalg.eigh(cov)
        XT_y = safe_sparse_dot(X.T, y, dense_output=True)
        XT_sqrt_sw = safe_sparse_dot(X.T, sqrt_sw, dense_output=True)
        if self.fit_intercept and sparse.issparse(X):
            # centered matrix = X - sqrt_sw X_mean'
            if y.ndim == 2:
                XT_y -= X_mean[:, None] * (sqrt_sw @ y)
            else:
                XT_y -= X_mean * (sqrt_sw @ y)
            XT_sqrt_sw -= X_mean * (sqrt_sw @ sqrt_sw)
        return eigvals, V, X, X_mean, XT_y, XT_sqrt_sw