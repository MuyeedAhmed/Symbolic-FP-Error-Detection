    def _eigen_decompose_gram(self, X, X_mean, y, sqrt_sw):
        """Eigendecomposition of Gram matrix X X'"""
        xp, is_array_api = get_namespace(X)
        K = self._compute_gram(X, X_mean, sqrt_sw)
        eigvals, Q = xp.linalg.eigh(K)
        QT_y = Q.T @ y
        QT_sqrt_sw = Q.T @ sqrt_sw
        XT = X.T
        return eigvals, Q, QT_y, QT_sqrt_sw, XT, X_mean