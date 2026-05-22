def test_decision_path(name):
    X = iris.data
    y = iris.target
    n_samples = X.shape[0]

    TreeEstimator = ALL_TREES[name]
    est = TreeEstimator(random_state=0, max_depth=2)
    est.fit(X, y)

    node_indicator_csr = est.decision_path(X)
    node_indicator = node_indicator_csr.toarray()
    assert node_indicator.shape == (n_samples, est.tree_.node_count)

    # Assert that leaves index are correct
    leaves = est.apply(X)
    leave_indicator = [node_indicator[i, j] for i, j in enumerate(leaves)]
    assert_array_almost_equal(leave_indicator, np.ones(shape=n_samples))

    # Ensure only one leave node per sample
    all_leaves = est.tree_.children_left == TREE_LEAF
    assert_array_almost_equal(
        np.dot(node_indicator, all_leaves), np.ones(shape=n_samples)
    )

    # Ensure max depth is consistent with sum of indicator
    max_depth = node_indicator.sum(axis=1).max()
    assert est.tree_.max_depth <= max_depth