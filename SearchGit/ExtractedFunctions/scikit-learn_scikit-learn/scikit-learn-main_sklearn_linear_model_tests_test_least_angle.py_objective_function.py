    def objective_function(coef):
        return 1.0 / (2.0 * len(X)) * linalg.norm(
            y - np.dot(X, coef)
        ) ** 2 + alpha * linalg.norm(coef, 1)