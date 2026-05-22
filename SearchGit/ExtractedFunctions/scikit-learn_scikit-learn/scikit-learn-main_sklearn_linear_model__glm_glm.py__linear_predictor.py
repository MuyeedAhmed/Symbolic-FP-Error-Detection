    def _linear_predictor(self, X):
        """Compute the linear_predictor = `X @ coef_ + intercept_`.

        Note that we often use the term raw_prediction instead of linear predictor.

        Parameters
        ----------
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Samples.

        Returns
        -------
        y_pred : array of shape (n_samples,)
            Returns predicted values of linear predictor.
        """
        xp, _ = get_namespace(X)
        check_is_fitted(self)
        X = validate_data(
            self,
            X,
            accept_sparse=["csr", "csc", "coo"],
            dtype=[xp.float64, xp.float32],
            ensure_2d=True,
            allow_nd=False,
            reset=False,
        )
        return X @ self.coef_ + self.intercept_