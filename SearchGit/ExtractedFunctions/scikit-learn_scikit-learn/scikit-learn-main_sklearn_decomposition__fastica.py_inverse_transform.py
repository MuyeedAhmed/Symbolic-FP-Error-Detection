    def inverse_transform(self, X, copy=True):
        """Transform the sources back to the mixed data (apply mixing matrix).

        Parameters
        ----------
        X : array-like of shape (n_samples, n_components)
            Sources, where `n_samples` is the number of samples
            and `n_components` is the number of components.
        copy : bool, default=True
            If False, data passed to fit are overwritten. Defaults to True.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Reconstructed data obtained with the mixing matrix.
        """
        check_is_fitted(self)

        X = check_array(X, copy=(copy and self.whiten), dtype=[np.float64, np.float32])
        X = np.dot(X, self.mixing_.T)
        if self.whiten:
            X += self.mean_

        return X