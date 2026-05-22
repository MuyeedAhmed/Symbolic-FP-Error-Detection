    def fit(self, X, y=None):
        """Fit estimator to data.

        Samples a subset of training points, computes kernel
        on these and computes normalization matrix.

        Parameters
        ----------
        X : array-like, shape (n_samples, n_features)
            Training data, where `n_samples` is the number of samples
            and `n_features` is the number of features.

        y : array-like, shape (n_samples,) or (n_samples, n_outputs), \
                default=None
            Target values (None for unsupervised transformations).

        Returns
        -------
        self : object
            Returns the instance itself.
        """
        xp, _, device = get_namespace_and_device(X)
        X = validate_data(self, X, accept_sparse="csr")
        rnd = check_random_state(self.random_state)
        n_samples = X.shape[0]

        # get basis vectors
        if self.n_components > n_samples:
            # XXX should we just bail?
            n_components = n_samples
            warnings.warn(
                "n_components > n_samples. This is not possible.\n"
                "n_components was set to n_samples, which results"
                " in inefficient evaluation of the full kernel."
            )

        else:
            n_components = self.n_components
        n_components = min(n_samples, n_components)
        inds = rnd.permutation(n_samples)
        basis_inds = xp.asarray(inds[:n_components], dtype=xp.int64, device=device)
        if sp.issparse(X):
            basis = X[basis_inds]
        else:
            basis = _safe_indexing(X, basis_inds, axis=0)

        basis_kernel = pairwise_kernels(
            basis,
            metric=self.kernel,
            filter_params=True,
            n_jobs=self.n_jobs,
            **self._get_kernel_params(),
        )

        # sqrt of kernel matrix on basis vectors
        _, _, dtype = _find_floating_dtype_allow_sparse(basis_kernel, Y=None, xp=xp)
        basis_kernel = xp.asarray(basis_kernel, dtype=dtype, device=device)
        U, S, V = xp.linalg.svd(basis_kernel)
        S = xp.clip(S, 1e-12, None)
        self.normalization_ = U / xp.sqrt(S) @ V
        self.components_ = basis
        self.component_indices_ = basis_inds
        self._n_features_out = n_components
        return self