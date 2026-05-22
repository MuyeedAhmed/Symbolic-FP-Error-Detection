    def fit(self, X, y=None):
        """Fit a Minimum Covariance Determinant with the FastMCD algorithm.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : Ignored
            Not used, present for API consistency by convention.

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        X = validate_data(self, X, ensure_min_samples=2, estimator="MinCovDet")
        random_state = check_random_state(self.random_state)
        n_samples, n_features = X.shape
        # check that the empirical covariance is full rank
        if (linalg.svdvals(np.dot(X.T, X)) > 1e-8).sum() != n_features:
            warnings.warn(
                "The covariance matrix associated to your dataset is not full rank"
            )
        # compute and store raw estimates
        raw_location, raw_covariance, raw_support, raw_dist = fast_mcd(
            X,
            support_fraction=self.support_fraction,
            cov_computation_method=self._nonrobust_covariance,
            random_state=random_state,
        )
        if self.assume_centered:
            raw_location = np.zeros(n_features)
            raw_covariance = self._nonrobust_covariance(
                X[raw_support], assume_centered=True
            )
            # get precision matrix in an optimized way
            precision = linalg.pinvh(raw_covariance)
            raw_dist = np.sum(np.dot(X, precision) * X, 1)
        self.raw_location_ = raw_location
        self.raw_covariance_ = raw_covariance
        self.raw_support_ = raw_support
        self.location_ = raw_location
        self.support_ = raw_support
        self.dist_ = raw_dist
        # obtain consistency at normal models
        self.correct_covariance(X)
        # re-weight estimator
        self.reweight_covariance(X)

        return self