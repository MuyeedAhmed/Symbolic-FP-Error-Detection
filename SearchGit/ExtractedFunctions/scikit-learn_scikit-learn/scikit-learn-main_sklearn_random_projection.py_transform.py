    def transform(self, X):
        """Project the data by using matrix product with the random matrix.

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_features)
            The input data to project into a smaller dimensional space.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            Projected array.
        """
        check_is_fitted(self)
        X = validate_data(
            self,
            X,
            accept_sparse=["csr", "csc"],
            reset=False,
            dtype=[np.float64, np.float32],
        )

        return _align_api_if_sparse(X @ self.components_.T)