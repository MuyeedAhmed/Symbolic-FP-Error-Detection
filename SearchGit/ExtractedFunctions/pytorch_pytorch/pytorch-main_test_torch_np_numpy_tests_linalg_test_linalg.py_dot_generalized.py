def dot_generalized(a, b):
    a = asarray(a)
    if a.ndim >= 3:
        if a.ndim == b.ndim:
            # matrix x matrix
            new_shape = a.shape[:-1] + b.shape[-1:]
        elif a.ndim == b.ndim + 1:
            # matrix x vector
            new_shape = a.shape[:-1]
        else:
            raise ValueError("Not implemented...")
        r = np.empty(new_shape, dtype=np.common_type(a, b))
        for c in itertools.product(*map(range, a.shape[:-2])):
            r[c] = dot(a[c], b[c])
        return r
    else:
        return dot(a, b)