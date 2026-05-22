def test_underflow_or_overflow():
    with np.errstate(all="raise"):
        # Generate some weird data with hugely unscaled features
        rng = np.random.RandomState(0)
        n_samples = 100
        n_features = 10

        X = rng.normal(size=(n_samples, n_features))
        X[:, :2] *= 1e300
        assert np.isfinite(X).all()

        # Use MinMaxScaler to scale the data without introducing a numerical
        # instability (computing the standard deviation naively is not possible
        # on this data)
        X_scaled = MinMaxScaler().fit_transform(X)
        assert np.isfinite(X_scaled).all()

        # Define a ground truth on the scaled data
        ground_truth = rng.normal(size=n_features)
        y = (np.dot(X_scaled, ground_truth) > 0.0).astype(np.int32)
        assert_array_equal(np.unique(y), [0, 1])

        model = SGDClassifier(alpha=0.1, loss="squared_hinge", max_iter=500)

        # smoke test: model is stable on scaled data
        model.fit(X_scaled, y)
        assert np.isfinite(model.coef_).all()

        # model is numerically unstable on unscaled data
        msg_regxp = (
            r"Floating-point under-/overflow occurred at epoch #.*"
            " Scaling input data with StandardScaler or MinMaxScaler"
            " might help."
        )
        with pytest.raises(ValueError, match=msg_regxp):
            model.fit(X, y)