def make_fullrank_matrices_with_distinct_singular_values(*shape, device, dtype, requires_grad=False):
    with torch.no_grad():
        t = make_tensor(shape, device=device, dtype=dtype)
        u, _, vh = torch.linalg.svd(t, full_matrices=False)
        real_dtype = t.real.dtype if t.dtype.is_complex else t.dtype
        k = min(shape[-1], shape[-2])
        # We choose the singular values to be "around one"
        # This is to make the matrix well conditioned
        # s = [2, 3, ..., k+1]
        s = torch.arange(2, k + 2, dtype=real_dtype, device=device)
        # s = [2, -3, 4, ..., (-1)^k k+1]
        s[1::2] *= -1.
        # 1 + 1/s so that the singular values are in the range [2/3, 3/2]
        # This gives a condition number of 9/4, which should be good enough
        s.reciprocal_().add_(1.)
        # Note that the singular values need not be ordered in an SVD so
        # we don't need to sort S
        x = (u * s.to(u.dtype)) @ vh
    x.requires_grad_(requires_grad)
    return x