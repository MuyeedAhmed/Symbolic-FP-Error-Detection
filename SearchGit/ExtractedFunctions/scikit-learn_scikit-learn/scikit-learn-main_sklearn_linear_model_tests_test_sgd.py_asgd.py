def asgd(klass, X, y, eta, alpha, weight_init=None, intercept_init=0.0):
    if weight_init is None:
        weights = np.zeros(X.shape[1])
    else:
        weights = weight_init

    average_weights = np.zeros(X.shape[1])
    intercept = intercept_init
    average_intercept = 0.0
    decay = 1.0

    # sparse data has a fixed decay of .01
    if klass in (SparseSGDClassifier, SparseSGDRegressor):
        decay = 0.01

    for i, entry in enumerate(X):
        p = np.dot(entry, weights)
        p += intercept
        gradient = p - y[i]
        weights *= 1.0 - (eta * alpha)
        weights += -(eta * gradient * entry)
        intercept += -(eta * gradient) * decay

        average_weights *= i
        average_weights += weights
        average_weights /= i + 1.0

        average_intercept *= i
        average_intercept += intercept
        average_intercept /= i + 1.0

    return average_weights, average_intercept