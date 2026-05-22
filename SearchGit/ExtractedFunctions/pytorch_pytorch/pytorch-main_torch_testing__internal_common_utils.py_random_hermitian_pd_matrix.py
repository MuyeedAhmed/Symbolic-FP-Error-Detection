def random_hermitian_pd_matrix(matrix_size, *batch_dims, dtype, device):
    """
    Returns a batch of random Hermitian positive-definite matrices.
    The shape of the result is batch_dims + (matrix_size, matrix_size)
    The following example creates a tensor of size 2 x 4 x 3 x 3
    >>> # xdoctest: +SKIP("undefined variables")
    >>> matrices = random_hermitian_pd_matrix(3, 2, 4, dtype=dtype, device=device)
    """
    A = torch.randn(*(batch_dims + (matrix_size, matrix_size)),
                    dtype=dtype, device=device)
    return A @ A.mH + torch.eye(matrix_size, dtype=dtype, device=device)