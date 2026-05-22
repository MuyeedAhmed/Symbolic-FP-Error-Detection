    def transform(self, X):
        """
        Transform new points into embedding space.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training set.

        Returns
        -------
        X_new : ndarray of shape (n_samples, n_components)
            Returns the instance itself.

        Notes
        -----
        Because of scaling performed by this method, it is discouraged to use
        it together with methods that are not scale-invariant (like SVMs).
        """
        check_is_fitted(self)

        X = validate_data(self, X, reset=False)
        ind = self.nbrs_.kneighbors(
            X, n_neighbors=self.n_neighbors, return_distance=False
        )
        weights = barycenter_weights(X, self.nbrs_._fit_X, ind, reg=self.reg)
        X_new = np.empty((X.shape[0], self.n_components))
        for i in range(X.shape[0]):
            X_new[i] = np.dot(self.embedding_[ind[i]].T, weights[i])
        return X_new