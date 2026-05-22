    def inverse_transform(self, X):
        """Project data back to its original space.

        Returns an array X_original whose transform would be X. Note that even
        if X is sparse, X_original is dense: this may use a lot of RAM.

        If `compute_inverse_components` is False, the inverse of the components is
        computed during each call to `inverse_transform` which can be costly.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_components)
            Data to be transformed back.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Reconstructed data.
        """
        check_is_fitted(self)

        X = check_array(X, dtype=[np.float64, np.float32], accept_sparse=("csr", "csc"))

        if self.compute_inverse_components:
            return _align_api_if_sparse(X @ self.inverse_components_.T)

        inverse_components = self._compute_inverse_components()
        return _align_api_if_sparse(X @ inverse_components.T)