def _laplace(m, d):
    return lambda v: v * d[:, np.newaxis] - m @ v