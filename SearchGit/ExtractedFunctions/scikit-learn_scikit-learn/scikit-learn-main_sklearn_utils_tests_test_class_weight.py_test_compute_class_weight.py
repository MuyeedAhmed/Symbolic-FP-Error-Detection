def test_compute_class_weight():
    # Test (and demo) compute_class_weight.
    y = np.asarray([2, 2, 2, 3, 3, 4])
    classes = np.unique(y)

    cw = compute_class_weight("balanced", classes=classes, y=y)
    # total effect of samples is preserved
    class_counts = np.bincount(y)[2:]
    assert_almost_equal(np.dot(cw, class_counts), y.shape[0])
    assert cw[0] < cw[1] < cw[2]