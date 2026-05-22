    def _solve_svd(self, X):
        """SVD solver for Quadratic Discriminant Analysis.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.
        """
        n_samples, n_features = X.shape

        mean = X.mean(0)
        Xc = X - mean
        # Xc = U * S * V.T
        _, S, Vt = np.linalg.svd(Xc, full_matrices=False)
        scaling = (S**2) / (n_samples - 1)  # scalings are squared singular values
        scaling = ((1 - self.reg_param) * scaling) + self.reg_param
        rotation = Vt.T

        cov = None
        if self.store_covariance:
            # cov = V * (S^2 / (n-1)) * V.T
            cov = scaling * Vt.T @ Vt

        return scaling, rotation, cov