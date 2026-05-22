def test_kernel_pca_precomputed(global_random_seed):
    """Test that kPCA works with a precomputed kernel, for all solvers"""
    rng = np.random.RandomState(global_random_seed)
    X_fit = rng.random_sample((5, 4))
    X_pred = rng.random_sample((2, 4))

    for eigen_solver in ("dense", "arpack", "randomized"):
        X_kpca = (
            KernelPCA(4, eigen_solver=eigen_solver, random_state=0)
            .fit(X_fit)
            .transform(X_pred)
        )

        X_kpca2 = (
            KernelPCA(
                4, eigen_solver=eigen_solver, kernel="precomputed", random_state=0
            )
            .fit(np.dot(X_fit, X_fit.T))
            .transform(np.dot(X_pred, X_fit.T))
        )

        X_kpca_train = KernelPCA(
            4, eigen_solver=eigen_solver, kernel="precomputed", random_state=0
        ).fit_transform(np.dot(X_fit, X_fit.T))

        X_kpca_train2 = (
            KernelPCA(
                4, eigen_solver=eigen_solver, kernel="precomputed", random_state=0
            )
            .fit(np.dot(X_fit, X_fit.T))
            .transform(np.dot(X_fit, X_fit.T))
        )

        assert_array_almost_equal(np.abs(X_kpca), np.abs(X_kpca2))

        assert_array_almost_equal(np.abs(X_kpca_train), np.abs(X_kpca_train2))