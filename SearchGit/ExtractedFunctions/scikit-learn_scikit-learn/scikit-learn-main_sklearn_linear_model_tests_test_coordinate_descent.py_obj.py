    def obj(W, X, y, alpha):
        intercept = W[:, -1]
        W = W[:, :-1]
        l21_norm = np.sqrt(np.sum(W**2, axis=0)).sum()
        return (
            np.linalg.norm(Y - X @ W.T - intercept, ord="fro") ** 2 / (2 * n_samples)
            + alpha * l21_norm
        )