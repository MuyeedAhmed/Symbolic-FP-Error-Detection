def test_conjugate_gradient() -> None:
    """
    >>> test_conjugate_gradient()  # self running tests
    """
    # Create linear system with SPD matrix and known solution x_true.
    dimension = 3
    spd_matrix = _create_spd_matrix(dimension)
    rng = np.random.default_rng()
    x_true = rng.normal(size=(dimension, 1))
    b = np.dot(spd_matrix, x_true)

    # Numpy solution.
    x_numpy = np.linalg.solve(spd_matrix, b)

    # Our implementation.
    x_conjugate_gradient = conjugate_gradient(spd_matrix, b)

    # Ensure both solutions are close to x_true (and therefore one another).
    assert np.linalg.norm(x_numpy - x_true) <= 1e-6
    assert np.linalg.norm(x_conjugate_gradient - x_true) <= 1e-6