def test_cv_pipeline_precomputed():
    # Cross-validate a regression on four coplanar points with the same
    # value. Use precomputed kernel to ensure Pipeline with KernelCenterer
    # is treated as a pairwise operation.
    X = np.array([[3, 0, 0], [0, 3, 0], [0, 0, 3], [1, 1, 1]])
    y_true = np.ones((4,))
    K = X.dot(X.T)
    kcent = KernelCenterer()
    pipeline = Pipeline([("kernel_centerer", kcent), ("svr", SVR())])

    # did the pipeline set the pairwise attribute?
    assert pipeline.__sklearn_tags__().input_tags.pairwise

    # test cross-validation, score should be almost perfect
    # NB: this test is pretty vacuous -- it's mainly to test integration
    #     of Pipeline and KernelCenterer
    y_pred = cross_val_predict(pipeline, K, y_true, cv=2)
    assert_array_almost_equal(y_true, y_pred)