def test_compute_class_weight_balanced_sample_weight_equivalence():
    # Test with unbalanced and negative class labels for
    # equivalence between repeated and weighted samples

    classes = np.array([-2, -1, 0])
    y = np.asarray([-1, -1, 0, 0, -2, -2])
    sw = np.asarray([1, 0, 1, 1, 1, 2])

    y_rep = np.repeat(y, sw, axis=0)

    class_weights_weighted = compute_class_weight(
        "balanced", classes=classes, y=y, sample_weight=sw
    )
    class_weights_repeated = compute_class_weight("balanced", classes=classes, y=y_rep)
    assert len(class_weights_weighted) == len(classes)
    assert len(class_weights_repeated) == len(classes)

    class_counts_weighted = np.bincount(y + 2, weights=sw)
    class_counts_repeated = np.bincount(y_rep + 2)

    assert np.dot(class_weights_weighted, class_counts_weighted) == pytest.approx(
        np.dot(class_weights_repeated, class_counts_repeated)
    )

    assert_allclose(class_weights_weighted, class_weights_repeated)