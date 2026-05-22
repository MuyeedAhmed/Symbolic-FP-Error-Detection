def asgd_oneclass(klass, X, eta, nu, coef_init=None, offset_init=0.0):
    if coef_init is None:
        coef = np.zeros(X.shape[1])
    else:
        coef = coef_init

    average_coef = np.zeros(X.shape[1])
    offset = offset_init
    intercept = 1 - offset
    average_intercept = 0.0
    decay = 1.0

    # sparse data has a fixed decay of .01
    if klass == SparseSGDOneClassSVM:
        decay = 0.01

    for i, entry in enumerate(X):
        p = np.dot(entry, coef)
        p += intercept
        if p <= 1.0:
            gradient = -1
        else:
            gradient = 0
        coef *= max(0, 1.0 - eta * nu)
        coef += -(eta * gradient * entry)
        intercept += -(eta * (nu + gradient)) * decay

        average_coef *= i
        average_coef += coef
        average_coef /= i + 1.0

        average_intercept *= i
        average_intercept += intercept
        average_intercept /= i + 1.0

    return average_coef, 1 - average_intercept