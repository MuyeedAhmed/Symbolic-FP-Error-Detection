def gen_toy_problem_2d():
    random_state = np.random.RandomState(0)
    n_samples = 100
    # Linear model y = 5*x_1 + 10*x_2 + N(1, 0.1**2)
    X = random_state.normal(size=(n_samples, 2))
    w = np.array([5.0, 10.0])
    c = 1.0
    noise = 0.1 * random_state.normal(size=n_samples)
    y = np.dot(X, w) + c + noise
    # Add some outliers
    n_outliers = n_samples // 10
    ix = random_state.randint(0, n_samples, size=n_outliers)
    y[ix] = 50 * random_state.normal(size=n_outliers)
    return X, y, w, c