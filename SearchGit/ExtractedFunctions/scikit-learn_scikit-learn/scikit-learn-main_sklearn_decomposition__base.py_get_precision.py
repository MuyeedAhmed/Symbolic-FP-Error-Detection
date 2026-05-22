    def get_precision(self):
        """Compute data precision matrix with the generative model.

        Equals the inverse of the covariance but computed with
        the matrix inversion lemma for efficiency.

        Returns
        -------
        precision : array, shape=(n_features, n_features)
            Estimated precision of data.
        """
        xp, is_array_api_compliant = get_namespace(self.components_)

        n_features = self.components_.shape[1]

        # handle corner cases first
        if self.n_components_ == 0:
            return xp.eye(n_features) / self.noise_variance_

        if is_array_api_compliant:
            linalg_inv = xp.linalg.inv
        else:
            linalg_inv = linalg.inv

        if self.noise_variance_ == 0.0:
            return linalg_inv(self.get_covariance())

        # Get precision using matrix inversion lemma
        components_ = self.components_
        exp_var = self.explained_variance_
        if self.whiten:
            components_ = components_ * xp.sqrt(exp_var[:, np.newaxis])
        exp_var_diff = exp_var - self.noise_variance_
        exp_var_diff = xp.where(
            exp_var > self.noise_variance_,
            exp_var_diff,
            xp.asarray(0.0, device=device(exp_var)),
        )
        precision = components_ @ components_.T / self.noise_variance_
        _add_to_diagonal(precision, 1.0 / exp_var_diff, xp)
        precision = components_.T @ linalg_inv(precision) @ components_
        precision /= -(self.noise_variance_**2)
        _add_to_diagonal(precision, 1.0 / self.noise_variance_, xp)
        return precision