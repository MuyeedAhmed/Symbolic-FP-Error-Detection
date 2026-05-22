def sample_inputs_linalg_cholesky_inverse(
    op_info, device, dtype, requires_grad=False, **kwargs
):
    from torch.testing._internal.common_utils import random_well_conditioned_matrix

    # Cholesky factorization is for positive-definite matrices
    single_well_conditioned_matrix = random_well_conditioned_matrix(
        S, S, dtype=dtype, device=device
    )
    batch_well_conditioned_matrices = random_well_conditioned_matrix(
        2, S, S, dtype=dtype, device=device
    )
    single_pd = single_well_conditioned_matrix @ single_well_conditioned_matrix.mH
    batch_pd = batch_well_conditioned_matrices @ batch_well_conditioned_matrices.mH

    inputs = (
        torch.zeros(0, 0, dtype=dtype, device=device),  # 0x0 matrix
        torch.zeros(0, 2, 2, dtype=dtype, device=device),  # zero batch of matrices
        single_pd,
        batch_pd,
    )
    test_cases = (torch.linalg.cholesky(a, upper=False) for a in inputs)
    for l in test_cases:
        # generated lower-triangular samples
        l.requires_grad = requires_grad
        yield SampleInput(l)  # upper=False by default
        yield SampleInput(
            l.detach().clone().requires_grad_(requires_grad), kwargs=dict(upper=False)
        )

        # generate upper-triangular inputs
        u = l.detach().clone().mT.contiguous().requires_grad_(requires_grad)
        yield SampleInput(u, kwargs=dict(upper=True))