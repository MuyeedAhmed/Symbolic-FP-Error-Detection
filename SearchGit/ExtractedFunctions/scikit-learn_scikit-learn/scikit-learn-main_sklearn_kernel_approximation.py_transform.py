    def transform(self, X):
        """Apply feature map to X.

        Computes an approximate feature map using the kernel
        between some training points and X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features)
            Data to transform.

        Returns
        -------
        X_transformed : ndarray of shape (n_samples, n_components)
            Transformed data.
        """
        check_is_fitted(self)

        xp, _, device = get_namespace_and_device(X)
        X = validate_data(self, X, accept_sparse="csr", reset=False)

        kernel_params = self._get_kernel_params()
        embedded = pairwise_kernels(
            X,
            self.components_,
            metric=self.kernel,
            filter_params=True,
            n_jobs=self.n_jobs,
            **kernel_params,
        )
        dtype = _find_matching_floating_dtype(embedded, xp=xp)
        embedded = xp.asarray(embedded, dtype=dtype, device=device)
        return embedded @ self.normalization_.T