    def _decision_function(self, X):
        check_is_fitted(self)

        X = validate_data(self, X, accept_sparse=["csr", "csc", "coo"], reset=False)
        coef_ = self.coef_
        if coef_.ndim == 1:
            return X @ coef_ + self.intercept_
        else:
            return X @ coef_.T + self.intercept_