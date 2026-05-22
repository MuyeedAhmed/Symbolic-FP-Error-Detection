def _create_spd_matrix(dimension: int) -> Any:
    """
    Returns a symmetric positive definite matrix given a dimension.

    Input:
    dimension gives the square matrix dimension.

    Output:
    spd_matrix is an diminesion x dimensions symmetric positive definite (SPD) matrix.

    >>> import numpy as np
    >>> dimension = 3
    >>> spd_matrix = _create_spd_matrix(dimension)
    >>> _is_matrix_spd(spd_matrix)
    True
    """
    rng = np.random.default_rng()
    random_matrix = rng.normal(size=(dimension, dimension))
    spd_matrix = np.dot(random_matrix, random_matrix.T)
    assert _is_matrix_spd(spd_matrix)
    return spd_matrix