def test_compute_class_weight_balanced_unordered():
    # Test compute_class_weight when classes are unordered
    classes = np.array([1, 0, 3])
    y = np.asarray([1, 0, 0, 3, 3, 3])

    cw = compute_class_weight("balanced", classes=classes, y=y)
    class_counts = np.bincount(y)[classes]
    assert_almost_equal(np.dot(cw, class_counts), y.shape[0])
    assert_array_almost_equal(cw, [2.0, 1.0, 2.0 / 3])