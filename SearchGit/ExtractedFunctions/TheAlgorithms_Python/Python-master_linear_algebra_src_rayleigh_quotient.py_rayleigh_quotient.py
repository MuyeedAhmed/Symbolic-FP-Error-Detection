def rayleigh_quotient(a: np.ndarray, v: np.ndarray) -> Any:
    """
    Returns the Rayleigh quotient of a Hermitian matrix A and
    vector v.
    >>> import numpy as np
    >>> A = np.array([
    ... [1,  2, 4],
    ... [2,  3,  -1],
    ... [4, -1,  1]
    ... ])
    >>> v = np.array([
    ... [1],
    ... [2],
    ... [3]
    ... ])
    >>> rayleigh_quotient(A, v)
    array([[3.]])
    """
    v_star = v.conjugate().T
    v_star_dot = v_star.dot(a)
    assert isinstance(v_star_dot, np.ndarray)
    return (v_star_dot.dot(v)) / (v_star.dot(v))