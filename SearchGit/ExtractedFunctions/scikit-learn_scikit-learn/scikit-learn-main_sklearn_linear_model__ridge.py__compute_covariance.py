    def _compute_covariance(self, X, X_mean, sqrt_sw):
        """Computes covariance matrix X' X with possible centering.

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
        covariance : ndarray of shape (n_features, n_features)
            The covariance matrix.

        Notes
        -----
        When self.fit_intercept is False no centering is done.

        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute X' X.

        When X is sparse it has not been centered in preprocessing, but
        it has been scaled by sqrt_sw. The centered X is never actually
        computed because centering would break the sparsity of X.
        """
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            return safe_sparse_dot(X.T, X, dense_output=True)
        # X is sparse and fit_intercept is True
        # centered matrix = X - sqrt_sw X_mean'
        sw_sum = sqrt_sw @ sqrt_sw
        return (
            safe_sparse_dot(X.T, X, dense_output=True)
            - sw_sum * X_mean[:, None] * X_mean[None, :]
        )