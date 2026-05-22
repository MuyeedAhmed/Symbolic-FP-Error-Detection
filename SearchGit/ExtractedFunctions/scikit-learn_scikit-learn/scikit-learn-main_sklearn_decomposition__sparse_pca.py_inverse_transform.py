    def inverse_transform(self, X):
        """Transform data from the latent space to the original space.

        This inversion is an approximation due to the loss of information
        induced by the forward decomposition.

        .. versionadded:: 1.2

        Parameters
        ----------
        X : ndarray of shape (n_samples, n_components)
            Data in the latent space.

        Returns
        -------
        X_original : ndarray of shape (n_samples, n_features)
            Reconstructed data in the original space.
        """
        check_is_fitted(self)
        X = check_array(X)

        return (X @ self.components_) + self.mean_