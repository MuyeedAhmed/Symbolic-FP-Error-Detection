    def transform(self, X):
        """Project data to maximize class separation.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Input data.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components) or \
            (n_samples, min(rank, n_components))
            Transformed data. In the case of the 'svd' solver, the shape
            is (n_samples, min(rank, n_components)).
        """
        if self.solver == "lsqr":
            raise NotImplementedError(
                "transform not implemented for 'lsqr' solver (use 'svd' or 'eigen')."
            )
        check_is_fitted(self)
        check_same_namespace(X, self, attribute="coef_", method="transform")
        X = validate_data(self, X, reset=False)

        if self.solver == "svd":
            X_new = (X - self.xbar_) @ self.scalings_
        elif self.solver == "eigen":
            X_new = X @ self.scalings_

        return X_new[:, : self._max_components]