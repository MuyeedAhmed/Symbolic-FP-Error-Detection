    def fit(self, X, y):
        """Fit model to data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training samples.

        y : array-like of shape (n_samples,) or (n_samples, n_targets)
            Targets.

        Returns
        -------
        self : object
            Fitted estimator.
        """
        check_consistent_length(X, y)
        X = validate_data(
            self,
            X,
            dtype=np.float64,
            force_writeable=True,
            copy=self.copy,
            ensure_min_samples=2,
        )
        y = check_array(
            y,
            input_name="y",
            dtype=np.float64,
            force_writeable=True,
            copy=self.copy,
            ensure_2d=False,
        )
        if y.ndim == 1:
            y = y.reshape(-1, 1)

        # we'll compute the SVD of the cross-covariance matrix = X.T.dot(y)
        # This matrix rank is at most min(n_samples, n_features, n_targets) so
        # n_components cannot be bigger than that.
        n_components = self.n_components
        rank_upper_bound = min(X.shape[0], X.shape[1], y.shape[1])
        if n_components > rank_upper_bound:
            raise ValueError(
                f"`n_components` upper bound is {rank_upper_bound}. "
                f"Got {n_components} instead. Reduce `n_components`."
            )

        X, y, self._x_mean, self._y_mean, self._x_std, self._y_std = _center_scale_xy(
            X, y, self.scale
        )

        # Compute SVD of cross-covariance matrix
        C = np.dot(X.T, y)
        U, s, Vt = svd(C, full_matrices=False)
        U = U[:, :n_components]
        Vt = Vt[:n_components]
        U, Vt = svd_flip(U, Vt)
        V = Vt.T

        self.x_weights_ = U
        self.y_weights_ = V
        self._n_features_out = self.x_weights_.shape[1]
        return self