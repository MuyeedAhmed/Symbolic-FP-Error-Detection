    def inverse_transform(self, X):
        """Transform data back to its original space.

        .. versionadded:: 0.18

        Parameters
        ----------
        X : {ndarray, sparse matrix} of shape (n_samples, n_components)
            Transformed data matrix.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Returns a data matrix of the original shape.
        """

        check_is_fitted(self)
        return X @ self.components_