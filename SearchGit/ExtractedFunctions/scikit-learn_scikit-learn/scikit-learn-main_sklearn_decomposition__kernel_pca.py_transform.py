    def transform(self, X):
        """Transform X.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training vector, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            Projection of X in the first principal components, where `n_samples`
            is the number of samples and `n_components` is the number of the components.
        """
        check_is_fitted(self)
        X = validate_data(self, X, accept_sparse="csr", reset=False)

        # Compute centered gram matrix between X and training data X_fit_
        K = self._centerer.transform(self._get_kernel(X, self.X_fit_))

        # scale eigenvectors (properly account for null-space for dot product)
        non_zeros = np.flatnonzero(self.eigenvalues_)
        scaled_alphas = np.zeros_like(self.eigenvectors_)
        scaled_alphas[:, non_zeros] = self.eigenvectors_[:, non_zeros] / np.sqrt(
            self.eigenvalues_[non_zeros]
        )

        # Project with a scalar product between K and the scaled eigenvectors
        return np.dot(K, scaled_alphas)