def test_x_none_gram_none_raises_value_error():
    # Test that lars_path with no X and Gram raises exception
    Xy = np.dot(X.T, y)
    with pytest.raises(ValueError, match="X and Gram cannot both be unspecified"):
        linear_model.lars_path(None, y, Gram=None, Xy=Xy)