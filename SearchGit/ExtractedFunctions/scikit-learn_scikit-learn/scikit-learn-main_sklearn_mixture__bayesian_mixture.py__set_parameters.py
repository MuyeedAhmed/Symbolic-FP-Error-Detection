    def _set_parameters(self, params, xp=None):
        (
            self.weight_concentration_,
            self.mean_precision_,
            self.means_,
            self.degrees_of_freedom_,
            self.covariances_,
            self.precisions_cholesky_,
        ) = params

        # Weights computation
        if self.weight_concentration_prior_type == "dirichlet_process":
            weight_dirichlet_sum = (
                self.weight_concentration_[0] + self.weight_concentration_[1]
            )
            tmp = self.weight_concentration_[1] / weight_dirichlet_sum
            self.weights_ = (
                self.weight_concentration_[0]
                / weight_dirichlet_sum
                * np.hstack((1, np.cumprod(tmp[:-1])))
            )
            self.weights_ /= np.sum(self.weights_)
        else:
            self.weights_ = self.weight_concentration_ / np.sum(
                self.weight_concentration_
            )

        # Precisions matrices computation
        if self.covariance_type == "full":
            self.precisions_ = np.array(
                [
                    np.dot(prec_chol, prec_chol.T)
                    for prec_chol in self.precisions_cholesky_
                ]
            )

        elif self.covariance_type == "tied":
            self.precisions_ = np.dot(
                self.precisions_cholesky_, self.precisions_cholesky_.T
            )
        else:
            self.precisions_ = self.precisions_cholesky_**2