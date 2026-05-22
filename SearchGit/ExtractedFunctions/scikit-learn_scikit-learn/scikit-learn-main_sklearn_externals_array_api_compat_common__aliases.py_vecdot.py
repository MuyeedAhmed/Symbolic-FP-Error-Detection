def vecdot(x1: Array, x2: Array, /, xp: Namespace, *, axis: int = -1) -> Array:
    if x1.shape[axis] != x2.shape[axis]:
        raise ValueError("x1 and x2 must have the same size along the given axis")

    if hasattr(xp, "broadcast_tensors"):
        _broadcast = xp.broadcast_tensors
    else:
        _broadcast = xp.broadcast_arrays

    x1_ = xp.moveaxis(x1, axis, -1)
    x2_ = xp.moveaxis(x2, axis, -1)
    x1_, x2_ = _broadcast(x1_, x2_)

    res = xp.conj(x1_[..., None, :]) @ x2_[..., None]
    return res[..., 0, 0]