    def transform(self, X, copy=True):
        """Recover the sources from X (apply the unmixing matrix).

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Data to transform, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        copy : bool, default=True
            If False, data passed to fit can be overwritten. Defaults to True.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            Estimated sources obtained by transforming the data with the
            estimated unmixing matrix.
        """
        check_is_fitted(self)

        X = validate_data(
            self,
            X,
            copy=(copy and self.whiten),
            dtype=[np.float64, np.float32],
            reset=False,
        )
        if self.whiten:
            X -= self.mean_

        return np.dot(X, self.components_.T)