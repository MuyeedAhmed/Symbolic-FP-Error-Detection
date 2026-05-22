def semi_sparse_linear(func, types, args=(), kwargs=None) -> torch.Tensor:
    if len(args) not in [2, 3]:
        raise AssertionError(f"expected 2 or 3 args, got {len(args)}")
    A, B = args[:2]
    bias = args[2] if len(args) == 3 else None

    shape = A.shape
    A_2d = A.view(-1, shape[-1])
    if bias is None:
        res = A_2d @ B.t()
    else:
        res = semi_sparse_addmm(
            func=None,
            types=None,
            args=[bias, A_2d, B.t()],
        )
    return res.view(*shape[:-1], -1)