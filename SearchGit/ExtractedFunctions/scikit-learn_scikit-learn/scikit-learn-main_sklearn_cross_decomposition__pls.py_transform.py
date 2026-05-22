    def transform(self, X, y=None):
        """
        Apply the dimensionality reduction.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Samples to be transformed.

        y : array-like of shape (n_samples,) or (n_samples, n_targets), \
                default=None
            Targets.

        Returns
        -------
        x_scores : array-like or tuple of array-like
            The transformed data `X_transformed` if `y is not None`,
            `(X_transformed, y_transformed)` otherwise.
        """
        check_is_fitted(self)
        X = validate_data(self, X, dtype=np.float64, reset=False)
        Xr = (X - self._x_mean) / self._x_std
        x_scores = np.dot(Xr, self.x_weights_)
        if y is not None:
            y = check_array(y, input_name="y", ensure_2d=False, dtype=np.float64)
            if y.ndim == 1:
                y = y.reshape(-1, 1)
            yr = (y - self._y_mean) / self._y_std
            y_scores = np.dot(yr, self.y_weights_)
            return x_scores, y_scores
        return x_scores