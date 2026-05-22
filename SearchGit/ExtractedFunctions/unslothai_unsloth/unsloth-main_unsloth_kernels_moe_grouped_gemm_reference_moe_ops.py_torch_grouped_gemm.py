def torch_grouped_gemm(X, W, m_sizes, transpose = True):
    """
    X: [M, K] if forward, else [M, N]
    W: [E, N, K]
    m_sizes: [E]

    Returns:
        Y: [M, N] if forward, else [M, K]
    """
    X = X.view(-1, X.shape[-1])
    M, K = X.shape

    assert m_sizes.ndim == 1
    E = m_sizes.shape[0]

    assert W.ndim == 3
    assert W.shape[0] == E

    N = W.shape[1]

    result = torch.zeros((M, N), dtype = X.dtype, device = X.device)

    m_start = 0
    for g in range(E):
        m_size = m_sizes[g]
        if m_size > 0:
            m_end = m_start + m_size

            # Extract group input
            # m_size x K
            X_g = X[m_start:m_end]
            # N x K
            W_g = W[g]

            # Y_g = X_g @ W_g.T -> [m_size, N]
            W_g = W_g.T if transpose else W_g
            Y_g = X_g @ W_g

            result[m_start:m_end] = Y_g

            m_start = m_end
    return result