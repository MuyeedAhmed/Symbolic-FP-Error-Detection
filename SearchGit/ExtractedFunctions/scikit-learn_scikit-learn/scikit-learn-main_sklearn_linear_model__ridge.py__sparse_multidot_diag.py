    def _sparse_multidot_diag(self, X, A, X_mean, sqrt_sw):
        """Compute the diagonal of X A X' with possible centering.

        Parameters
        ----------
        X : {ndarray, sparse matrix, sparse array} of shape (n_samples, n_features)
            The preprocessed design matrix.

        A : ndarray of shape (n_features, n_features)
            The inner matrix.

        X_mean : ndarray of shape (n_feature,)
            The weighted mean of X for each feature.

        sqrt_sw : ndarray of shape (n_samples,)
            Square roots of sample weights.

        Returns
        -------
        diag : np.ndarray, shape (n_samples,)
            The computed diagonal.

        Notes
        -----
        When self.fit_intercept is False no centering is done.

        When X is dense the centering has been done in preprocessing
        so the mean is 0 and we just compute diag(X A X').

        When X is sparse it has not been centered in preprocessing, but
        it has been scaled by sqrt_sw. The centered X is never actually
        computed because centering would break the sparsity of X.
        """
        xp, _ = get_namespace(X)
        XA = X @ A
        if sparse.isspmatrix(X):
            # sparse matrix use multiply for element wise multiplication
            XAX = np.ravel(X.multiply(XA).sum(axis=1))
        else:
            XAX = xp.sum(XA * X, axis=1)
        center = self.fit_intercept and sparse.issparse(X)
        if not center:
            # in this case centering has been done in preprocessing
            # or we are not fitting an intercept.
            return XAX
        # X is sparse and fit_intercept is True
        # centered matrix = X - sqrt_sw X_mean'
        XA_Xm = XA @ X_mean
        A_Xm = A @ X_mean
        sw = sqrt_sw * sqrt_sw
        return XAX - 2 * sqrt_sw * XA_Xm + sw * (X_mean @ A_Xm)