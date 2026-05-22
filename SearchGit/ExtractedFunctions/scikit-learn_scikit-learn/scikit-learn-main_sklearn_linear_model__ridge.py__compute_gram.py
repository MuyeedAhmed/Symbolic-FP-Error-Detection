    def _compute_gram(self, X, X_mean, sqrt_sw):
        """Computes the Gram matrix X X' with possible centering.

        Parameters
        ----------
        X : {ndarray, sparse matrix, sparse array} of shape (n_samples, n_features)
            The preprocessed design matrix.

        X_mean : ndarray of shape (n_feature,)
            The weighted mean of X for each feature.

        sqrt_sw : ndarray of shape (n_samples,)
            Square roots of sample weights.

        Returns
        -------
        gram : ndarray of shape (n_samples, n_samples)
            The Gram matrix.

        Notes
        -----
        When self.fit_intercept is False no centering is done.

        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute X X'.

        When X is sparse it has not been centered in preprocessing, but
        it has been scaled by sqrt_sw. The centered X is never actually
        computed because centering would break the sparsity of X.
        """
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            return safe_sparse_dot(X, X.T, dense_output=True)
        # X is sparse and fit_intercept is True
        # centered matrix = X - sqrt_sw X_mean'
        X_Xm = safe_sparse_dot(X, X_mean, dense_output=True)
        return (
            safe_sparse_dot(X, X.T, dense_output=True)
            - X_Xm[:, None] * sqrt_sw[None, :]
            - sqrt_sw[:, None] * X_Xm[None, :]
            + (X_mean @ X_mean) * sqrt_sw[:, None] * sqrt_sw[None, :]
        )