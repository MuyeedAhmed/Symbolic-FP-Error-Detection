    def inverse_transform(self, X):
        """Transform data back to its original space.

        In other words, return an input `X_original` whose transform would be X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_components)
            New data, where `n_samples` is the number of samples
            and `n_components` is the number of components.

        Returns
        -------
        X_original : array-like of shape (n_samples, n_features)
            Original data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        Notes
        -----
        If whitening is enabled, inverse_transform will compute the
        exact inverse operation, which includes reversing whitening.
        """
        xp, _ = get_namespace(X, self.components_, self.explained_variance_)

        check_is_fitted(self)

        X = check_array(X, input_name="X", dtype=[xp.float64, xp.float32])

        if self.whiten:
            scaled_components = (
                xp.sqrt(self.explained_variance_[:, np.newaxis]) * self.components_
            )
            return X @ scaled_components + self.mean_
        else:
            return X @ self.components_ + self.mean_