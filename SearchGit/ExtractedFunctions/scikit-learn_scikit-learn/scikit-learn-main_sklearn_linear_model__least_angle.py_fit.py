    def fit(self, X, y, copy_X=None):
        """Fit the model using X, y as training data.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Training data.

        y : array-like of shape (n_samples,)
            Target values. Will be cast to X's dtype if necessary.

        copy_X : bool, default=None
            If provided, this parameter will override the choice
            of copy_X made at instance creation.
            If ``True``, X will be copied; else, it may be overwritten.

        Returns
        -------
        self : object
            Returns an instance of self.
        """
        if copy_X is None:
            copy_X = self.copy_X
        X, y = validate_data(self, X, y, force_writeable=True, y_numeric=True)

        X, y, Xmean, ymean, _, _ = _preprocess_data(
            X, y, fit_intercept=self.fit_intercept, copy=copy_X
        )

        Gram = self.precompute

        alphas_, _, coef_path_, self.n_iter_ = lars_path(
            X,
            y,
            Gram=Gram,
            copy_X=copy_X,
            copy_Gram=True,
            alpha_min=0.0,
            method="lasso",
            verbose=self.verbose,
            max_iter=self.max_iter,
            eps=self.eps,
            return_n_iter=True,
            positive=self.positive,
        )

        n_samples = X.shape[0]

        if self.criterion == "aic":
            criterion_factor = 2
        elif self.criterion == "bic":
            criterion_factor = log(n_samples)
        else:
            raise ValueError(
                f"criterion should be either bic or aic, got {self.criterion!r}"
            )

        residuals = y[:, np.newaxis] - np.dot(X, coef_path_)
        residuals_sum_squares = np.sum(residuals**2, axis=0)
        degrees_of_freedom = np.zeros(coef_path_.shape[1], dtype=int)
        for k, coef in enumerate(coef_path_.T):
            mask = np.abs(coef) > np.finfo(coef.dtype).eps
            if not np.any(mask):
                continue
            # get the number of degrees of freedom equal to:
            # Xc = X[:, mask]
            # Trace(Xc * inv(Xc.T, Xc) * Xc.T) ie the number of non-zero coefs
            degrees_of_freedom[k] = np.sum(mask)

        self.alphas_ = alphas_

        if self.noise_variance is None:
            self.noise_variance_ = self._estimate_noise_variance(
                X, y, positive=self.positive
            )
        else:
            self.noise_variance_ = self.noise_variance

        self.criterion_ = (
            n_samples * np.log(2 * np.pi * self.noise_variance_)
            + residuals_sum_squares / self.noise_variance_
            + criterion_factor * degrees_of_freedom
        )
        n_best = np.argmin(self.criterion_)

        self.alpha_ = alphas_[n_best]
        self.coef_ = coef_path_[:, n_best]
        self._set_intercept(Xmean, ymean)
        return self