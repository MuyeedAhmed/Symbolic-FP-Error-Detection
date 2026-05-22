def _inv(x):
    x = backend.convert_to_tensor(x)
    _assert_2d(x)
    _assert_square(x)
    return backend.linalg.inv(x)