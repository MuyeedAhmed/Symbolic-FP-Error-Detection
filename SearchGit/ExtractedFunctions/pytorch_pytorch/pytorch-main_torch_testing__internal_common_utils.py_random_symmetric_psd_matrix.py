def random_symmetric_psd_matrix(l, *batches, **kwargs):
    """
    Returns a batch of random symmetric positive-semi-definite matrices.
    The shape of the result is batch_dims + (matrix_size, matrix_size)
    The following example creates a tensor of size 2 x 4 x 3 x 3
    >>> # xdoctest: +SKIP("undefined variables")
    >>> matrices = random_symmetric_psd_matrix(3, 2, 4, dtype=dtype, device=device)
    """
    dtype = kwargs.get('dtype', torch.double)
    device = kwargs.get('device', 'cpu')
    A = torch.randn(*(batches + (l, l)), dtype=dtype, device=device)
    return A @ A.mT