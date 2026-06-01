        def update_coeff(X, y, coef_, alpha_, keep_lambda, sigma_):
            coef_[keep_lambda] = alpha_ * np.linalg.multi_dot(
                [sigma_, X[:, keep_lambda].T, y]
            )
            return coef_