    def predict(self, X, copy=True):
        """Predict targets of given samples.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Samples.

        copy : bool, default=True
            Whether to copy `X` or perform in-place normalization.

        Returns
        -------
        y_pred : ndarray of shape (n_samples,) or (n_samples, n_targets)
            Returns predicted values.

        Notes
        -----
        This call requires the estimation of a matrix of shape
        `(n_features, n_targets)`, which may be an issue in high dimensional
        space.
        """
        check_is_fitted(self)
        X = validate_data(self, X, copy=copy, dtype=FLOAT_DTYPES, reset=False)
        # Only center X but do not scale it since the coefficients are already scaled
        X -= self._x_mean
        y_pred = X @ self.coef_.T + self.intercept_
        return y_pred.ravel() if self._predict_1d else y_pred