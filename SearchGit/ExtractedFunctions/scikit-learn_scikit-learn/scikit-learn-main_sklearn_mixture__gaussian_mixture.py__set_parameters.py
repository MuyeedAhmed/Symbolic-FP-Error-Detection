    def _set_parameters(self, params, xp=None):
        xp, _, device_ = get_namespace_and_device(params, xp=xp)
        (
            self.weights_,
            self.means_,
            self.covariances_,
            self.precisions_cholesky_,
        ) = params

        # Attributes computation
        if self.covariance_type == "full":
            self.precisions_ = xp.empty_like(self.precisions_cholesky_, device=device_)
            for k in range(self.precisions_cholesky_.shape[0]):
                prec_chol = self.precisions_cholesky_[k, :, :]
                self.precisions_[k, :, :] = prec_chol @ prec_chol.T

        elif self.covariance_type == "tied":
            self.precisions_ = self.precisions_cholesky_ @ self.precisions_cholesky_.T

        else:
            self.precisions_ = self.precisions_cholesky_**2