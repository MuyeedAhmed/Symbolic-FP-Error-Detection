def random_well_conditioned_matrix(*shape, dtype, device, mean=1.0, sigma=0.001):
    """
    Returns a random rectangular matrix (batch of matrices)
    with singular values sampled from a Gaussian with
    mean `mean` and standard deviation `sigma`.
    The smaller the `sigma`, the better conditioned
    the output matrix is.
    """
    primitive_dtype = {
        torch.float: torch.float,
        torch.double: torch.double,
        torch.cfloat: torch.float,
        torch.cdouble: torch.double
    }
    x = torch.rand(shape, dtype=dtype, device=device)
    m = x.size(-2)
    n = x.size(-1)
    u, _, vh = torch.linalg.svd(x, full_matrices=False)
    s = (torch.randn(*(shape[:-2] + (min(m, n),)), dtype=primitive_dtype[dtype], device=device) * sigma + mean) \
        .sort(-1, descending=True).values.to(dtype)
    return (u * s.unsqueeze(-2)) @ vh