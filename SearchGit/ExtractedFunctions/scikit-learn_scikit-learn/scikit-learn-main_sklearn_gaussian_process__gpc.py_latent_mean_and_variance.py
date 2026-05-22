    def latent_mean_and_variance(self, X):
        """Compute the mean and variance of the latent function values.

        Based on algorithm 3.2 of [RW2006]_, this function returns the latent
        mean (Line 4) and variance (Line 6) of the Gaussian process
        classification model.

        Note that this function is only supported for binary classification.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Query points where the GP is evaluated for classification.

        Returns
        -------
        latent_mean : array-like of shape (n_samples,)
            Mean of the latent function values at the query points.

        latent_var : array-like of shape (n_samples,)
            Variance of the latent function values at the query points.
        """
        check_is_fitted(self)

        # Based on Algorithm 3.2 of GPML
        K_star = self.kernel_(self.X_train_, X)  # K_star =k(x_star)
        latent_mean = K_star.T.dot(self.y_train_ - self.pi_)  # Line 4
        v = solve(self.L_, self.W_sr_[:, np.newaxis] * K_star)  # Line 5
        # Line 6 (compute np.diag(v.T.dot(v)) via einsum)
        latent_var = self.kernel_.diag(X) - np.einsum("ij,ij->j", v, v)

        return latent_mean, latent_var