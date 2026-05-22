def random_square_matrix_of_rank(l, rank, dtype=torch.double, device='cpu'):
    if rank > l:
        raise AssertionError(f"rank must be <= l, got rank={rank}, l={l}")
    A = torch.randn(l, l, dtype=dtype, device=device)
    u, s, vh = torch.linalg.svd(A, full_matrices=False)
    for i in range(l):
        if i >= rank:
            s[i] = 0
        elif s[i] == 0:
            s[i] = 1
    return (u * s.to(dtype).unsqueeze(-2)) @ vh