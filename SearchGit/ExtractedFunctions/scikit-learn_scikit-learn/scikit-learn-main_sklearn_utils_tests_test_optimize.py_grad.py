    def grad(x):
        return A.T @ (A @ x - y)