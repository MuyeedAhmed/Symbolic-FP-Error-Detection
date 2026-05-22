def random_hermitian_psd_matrix(matrix_size, *batch_dims, dtype=torch.double, device='cpu'):
    """
    Returns a batch of random Hermitian positive-semi-definite matrices.
    The shape of the result is batch_dims + (matrix_size, matrix_size)
    The following example creates a tensor of size 2 x 4 x 3 x 3
    >>> # xdoctest: +SKIP("undefined variables")
    >>> matrices = random_hermitian_psd_matrix(3, 2, 4, dtype=dtype, device=device)
    """
    A = torch.randn(*(batch_dims + (matrix_size, matrix_size)), dtype=dtype, device=device)
    return A @ A.mH