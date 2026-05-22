def random_matrix(rows, columns, *batch_dims, **kwargs):
    """Return rectangular matrix or batches of rectangular matrices.

    Parameters:
      dtype - the data type
      device - the device kind
      singular - when True, the output will be singular
    """
    dtype = kwargs.get('dtype', torch.double)
    device = kwargs.get('device', 'cpu')
    silent = kwargs.get("silent", False)
    singular = kwargs.get("singular", False)
    if silent and not torch._C.has_lapack:
        return torch.ones(rows, columns, dtype=dtype, device=device)

    A = torch.randn(batch_dims + (rows, columns), dtype=dtype, device=device)
    if A.numel() == 0:
        return A
    u, _, vh = torch.linalg.svd(A, full_matrices=False)
    k = min(rows, columns)
    s = torch.linspace(1 / (k + 1), 1, k, dtype=dtype, device=device)
    if singular:
        # make matrix singular
        s[k - 1] = 0
        if k > 2:
            # increase the order of singularity so that the pivoting
            # in LU factorization will be non-trivial
            s[0] = 0
    return (u * s.unsqueeze(-2)) @ vh