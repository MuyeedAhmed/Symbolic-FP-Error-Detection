    def ridge_loss(model, random_state=None, noise_scale=1e-8):
        intercept = model.intercept_
        if random_state is not None:
            rng = np.random.RandomState(random_state)
            coef = model.coef_ + rng.uniform(0, noise_scale, size=model.coef_.shape)
        else:
            coef = model.coef_

        return 0.5 * np.sum((y - X @ coef - intercept) ** 2) + 0.5 * alpha * np.sum(
            coef**2
        )