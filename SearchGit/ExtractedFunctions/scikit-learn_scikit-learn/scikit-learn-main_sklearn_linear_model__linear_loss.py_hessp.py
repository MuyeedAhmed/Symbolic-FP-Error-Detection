            def hessp(s):
                s = s.reshape((n_classes, -1), order="F")  # shape = (n_classes, n_dof)
                if self.fit_intercept:
                    s_intercept = s[:, -1]
                    s = s[:, :-1]  # shape = (n_classes, n_features)
                else:
                    s_intercept = 0
                tmp = X @ s.T + s_intercept  # X_{im} * s_k_m
                tmp -= (proba * tmp).sum(axis=1)[:, np.newaxis]  # - sum_l ..
                tmp *= proba  # * p_i_k
                if sample_weight is not None:
                    tmp *= sample_weight[:, np.newaxis]
                # hess_prod = empty_like(grad), but we ravel grad below and this
                # function is run after that.
                hess_prod = np.empty((n_classes, n_dof), dtype=weights.dtype, order="F")
                hess_prod[:, :n_features] = (tmp.T @ X) / sw_sum + l2_reg_strength * s
                if self.fit_intercept:
                    hess_prod[:, -1] = tmp.sum(axis=0) / sw_sum
                if coef.ndim == 1:
                    return hess_prod.ravel(order="F")
                else:
                    return hess_prod