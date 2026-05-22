    def transform(self, X):
        """Apply the learned transformation to the given data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Data samples.

        Returns
        -------
        X_embedded: ndarray of shape (n_samples, n_components)
            The data samples transformed.

        Raises
        ------
        NotFittedError
            If :meth:`fit` has not been called before.
        """

        check_is_fitted(self)
        X = validate_data(self, X, reset=False)

        return np.dot(X, self.components_.T)