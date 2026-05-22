def transform(velocity: float, event: np.ndarray | None = None) -> np.ndarray:
    """
    Calculate a Lorentz transformation for movement in the x direction given a
    velocity and a four-vector for an inertial reference frame

    If no four-vector is given, then calculate the transformation symbolically
    with variables
    >>> transform(29979245, np.array([1, 2, 3, 4]))
    array([ 3.01302757e+08, -3.01302729e+07,  3.00000000e+00,  4.00000000e+00])
    >>> transform(29979245)
    array([1.00503781498831*ct - 0.100503778816875*x,
           -0.100503778816875*ct + 1.00503781498831*x, 1.0*y, 1.0*z],
          dtype=object)
    >>> transform(19879210.2)
    array([1.0022057787097*ct - 0.066456172618675*x,
           -0.066456172618675*ct + 1.0022057787097*x, 1.0*y, 1.0*z],
          dtype=object)
    >>> transform(299792459, np.array([1, 1, 1, 1]))
    Traceback (most recent call last):
      ...
    ValueError: Speed must not exceed light speed 299,792,458 [m/s]!
    >>> transform(-1, np.array([1, 1, 1, 1]))
    Traceback (most recent call last):
      ...
    ValueError: Speed must be greater than or equal to 1!
    """
    # Ensure event is not empty
    if event is None:
        event = np.array([ct, x, y, z])  # Symbolic four vector
    else:
        event[0] *= c  # x0 is ct (speed of light * time)

    return transformation_matrix(velocity) @ event