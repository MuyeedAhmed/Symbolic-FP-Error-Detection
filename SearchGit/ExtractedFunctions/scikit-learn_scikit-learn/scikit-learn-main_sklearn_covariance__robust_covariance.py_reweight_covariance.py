    def reweight_covariance(self, data):
        """Re-weight raw Minimum Covariance Determinant estimates.

        Re-weight observations using Rousseeuw's method (equivalent to
        deleting outlying observations from the data set before
        computing location and covariance estimates) described
        in [RVDriessen]_.

        Corrects the re-weighted covariance to be consistent at the normal
        distribution, following [Croux1999]_.

        Parameters
        ----------
        data : array-like of shape (n_samples, n_features)
            The data matrix, with p features and n samples.
            The data set must be the one which was used to compute
            the raw estimates.

        Returns
        -------
        location_reweighted : ndarray of shape (n_features,)
            Re-weighted robust location estimate.

        covariance_reweighted : ndarray of shape (n_features, n_features)
            Re-weighted robust covariance estimate.

        support_reweighted : ndarray of shape (n_samples,), dtype=bool
            A mask of the observations that have been used to compute
            the re-weighted robust location and covariance estimates.

        References
        ----------

        .. [RVDriessen] A Fast Algorithm for the Minimum Covariance
            Determinant Estimator, 1999, American Statistical Association
            and the American Society for Quality, TECHNOMETRICS

        .. [Croux1999] Influence Function and Efficiency of the Minimum
            Covariance Determinant Scatter Matrix Estimator, 1999, Journal of
            Multivariate Analysis, Volume 71, Issue 2, Pages 161-190
        """
        n_samples, n_features = data.shape
        quantile_threshold = 0.025
        mask = self.dist_ < chi2(n_features).isf(quantile_threshold)
        if self.assume_centered:
            location_reweighted = np.zeros(n_features)
        else:
            location_reweighted = data[mask].mean(0)
        covariance_reweighted = self._nonrobust_covariance(
            data[mask], assume_centered=self.assume_centered
        )
        support_reweighted = np.zeros(n_samples, dtype=bool)
        support_reweighted[mask] = True
        # Parameter alpha as in [Croux1999] Eq. 4.2
        consistency_factor = _consistency_factor(
            n_features=n_features, alpha=1 - quantile_threshold
        )
        self._set_covariance(covariance_reweighted * consistency_factor)
        self.location_ = location_reweighted
        self.support_ = support_reweighted
        X_centered = data - self.location_
        self.dist_ = np.sum(np.dot(X_centered, self.get_precision()) * X_centered, 1)
        return location_reweighted, covariance_reweighted, support_reweighted