    def _svd_decompose_design_matrix(self, X, X_mean, y, sqrt_sw):
        """Reduced SVD decomposition of X"""
        xp, _ = get_namespace(X)
        # reduced svd
        U, singvals, VT = xp.linalg.svd(X, full_matrices=False)
        UT_y = U.T @ y
        UT_sqrt_sw = U.T @ sqrt_sw
        V = VT.T
        return singvals, U, V, UT_y, UT_sqrt_sw