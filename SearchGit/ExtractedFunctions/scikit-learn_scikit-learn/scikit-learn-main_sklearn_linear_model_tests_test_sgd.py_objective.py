    def objective(w, X, y, alpha):
        weights = w[:-1]
        intercept = w[-1]
        p = X @ weights + intercept
        z = p * y
        avg_loss = np.mean(np.maximum(hinge_threshold - z, 0.0))
        reg = 0.5 * alpha * weights @ weights
        obj = avg_loss + reg + intercept * alpha
        return obj