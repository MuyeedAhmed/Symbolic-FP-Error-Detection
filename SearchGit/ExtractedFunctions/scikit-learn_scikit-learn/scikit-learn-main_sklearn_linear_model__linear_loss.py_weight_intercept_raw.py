    def weight_intercept_raw(self, coef, X):
        """Helper function to get coefficients, intercept and raw_prediction.

        Parameters
        ----------
        coef : ndarray of shape (n_dof,), (n_classes, n_dof) or (n_classes * n_dof,)
            Coefficients of a linear model.
            If shape (n_classes * n_dof,), the classes of one feature are contiguous,
            i.e. one reconstructs the 2d-array via
            coef.reshape((n_classes, -1), order="F").
        X : {array-like, sparse matrix} of shape (n_samples, n_features)
            Training data.

        Returns
        -------
        weights : ndarray of shape (n_features,) or (n_classes, n_features)
            Coefficients without intercept term.
        intercept : float or ndarray of shape (n_classes,)
            Intercept terms.
        raw_prediction : ndarray of shape (n_samples,) or \
            (n_samples, n_classes)
        """
        weights, intercept = self.weight_intercept(coef)
        xp, _, device_ = get_namespace_and_device(X)

        # The `weights` and `intercept` are only converted internally to the
        # array API because the relevant `scipy.optimize` functions do not
        # currently support the array API and we have to ensure that the final
        # values returned to the respective `scipy.optimize` function are in
        # the `numpy` namespace.
        weights_xp = xp.asarray(weights, dtype=X.dtype, device=device_)
        intercept_xp = xp.asarray(intercept, dtype=X.dtype, device=device_)
        if not self.base_loss.is_multiclass:
            raw_prediction = X @ weights_xp + intercept_xp
        else:
            # weights has shape (n_classes, n_dof)
            raw_prediction = X @ weights_xp.T + intercept_xp

        return weights, intercept, raw_prediction