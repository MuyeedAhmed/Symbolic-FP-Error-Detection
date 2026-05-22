def inv(a: ArrayLike):
    a = _atleast_float_1(a)
    result = torch.linalg.inv(a)
    return result