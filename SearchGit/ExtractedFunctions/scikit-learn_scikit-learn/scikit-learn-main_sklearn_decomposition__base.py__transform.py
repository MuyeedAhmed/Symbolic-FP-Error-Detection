    def _transform(self, X, xp, x_is_centered=False):
        X_transformed = X @ self.components_.T
        if not x_is_centered:
            # Apply the centering after the projection.
            # For dense X this avoids copying or mutating the data passed by
            # the caller.
            # For sparse X it keeps sparsity and avoids having to wrap X into
            # a linear operator.
            X_transformed -= xp.reshape(self.mean_, (1, -1)) @ self.components_.T
        if self.whiten:
            # For some solvers (such as "arpack" and "covariance_eigh"), on
            # rank deficient data, some components can have a variance
            # arbitrarily close to zero, leading to non-finite results when
            # whitening. To avoid this problem we clip the variance below.
            scale = xp.sqrt(self.explained_variance_)
            min_scale = xp.finfo(scale.dtype).eps
            scale[scale < min_scale] = min_scale
            X_transformed /= scale
        return X_transformed