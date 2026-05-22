    def inverse_transform(self, X):
        """Transform X back to its original space.

        Returns an array X_original whose transform would be X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_components)
            New data.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Note that this is always a dense array.
        """
        X = check_array(X)
        return np.dot(X, self.components_)