    def _set_intercept(self, X_offset, y_offset, X_scale=None):
        """Set the intercept_"""
        xp, _ = get_namespace(X_offset, y_offset, X_scale)

        if self.fit_intercept:
            # We always want coef_.dtype=X.dtype. For instance, X.dtype can differ from
            # coef_.dtype if warm_start=True.
            self.coef_ = xp.astype(self.coef_, X_offset.dtype, copy=False)
            if X_scale is not None:
                self.coef_ = xp.divide(self.coef_, X_scale)

            if self.coef_.ndim == 1:
                self.intercept_ = y_offset - X_offset @ self.coef_
            else:
                self.intercept_ = y_offset - X_offset @ self.coef_.T

        else:
            self.intercept_ = 0.0