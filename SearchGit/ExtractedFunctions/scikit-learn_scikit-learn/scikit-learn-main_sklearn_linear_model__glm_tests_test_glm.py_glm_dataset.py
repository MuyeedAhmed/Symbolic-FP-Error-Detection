def glm_dataset(global_random_seed, request):
    """Dataset with GLM solutions, well conditioned X.

    This is inspired by ols_ridge_dataset in test_ridge.py.

    The construction is based on the SVD decomposition of X = U S V'.

    Parameters
    ----------
    type : {"long", "wide"}
        If "long", then n_samples > n_features.
        If "wide", then n_features > n_samples.
    model : a GLM model

    For "wide", we return the minimum norm solution:

        min ||w||_2 subject to w = argmin deviance(X, y, w)

    Note that the deviance is always minimized if y = inverse_link(X w) is possible to
    achieve, which it is in the wide data case. Therefore, we can construct the
    solution with minimum norm like (wide) OLS:

        min ||w||_2 subject to link(y) = raw_prediction = X w

    Returns
    -------
    model : GLM model
    X : ndarray
        Last column of 1, i.e. intercept.
    y : ndarray
    coef_unpenalized : ndarray
        Minimum norm solutions, i.e. min sum(loss(w)) (with minimum ||w||_2 in
        case of ambiguity)
        Last coefficient is intercept.
    coef_penalized : ndarray
        GLM solution with alpha=l2_reg_strength=1, i.e.
        min 1/n * sum(loss) + ||w[:-1]||_2^2.
        Last coefficient is intercept.
    l2_reg_strength : float
        Always equal 1.
    """
    data_type, model = request.param
    # Make larger dim more than double as big as the smaller one.
    # This helps when constructing singular matrices like (X, X).
    if data_type == "long":
        n_samples, n_features = 12, 4
    else:
        n_samples, n_features = 4, 12
    k = min(n_samples, n_features)
    rng = np.random.RandomState(global_random_seed)
    X = make_low_rank_matrix(
        n_samples=n_samples,
        n_features=n_features,
        effective_rank=k,
        tail_strength=0.1,
        random_state=rng,
    )
    X[:, -1] = 1  # last columns acts as intercept
    U, s, Vt = linalg.svd(X, full_matrices=False)
    assert np.all(s > 1e-3)  # to be sure
    assert np.max(s) / np.min(s) < 100  # condition number of X

    if data_type == "long":
        coef_unpenalized = rng.uniform(low=1, high=3, size=n_features)
        coef_unpenalized *= rng.choice([-1, 1], size=n_features)
        raw_prediction = X @ coef_unpenalized
    else:
        raw_prediction = rng.uniform(low=-3, high=3, size=n_samples)
        # minimum norm solution min ||w||_2 such that raw_prediction = X w:
        # w = X'(XX')^-1 raw_prediction = V s^-1 U' raw_prediction
        coef_unpenalized = Vt.T @ np.diag(1 / s) @ U.T @ raw_prediction

    linear_loss = LinearModelLoss(base_loss=model._get_loss(), fit_intercept=True)
    sw = np.full(shape=n_samples, fill_value=1 / n_samples)
    y = linear_loss.base_loss.link.inverse(raw_prediction)

    # Add penalty l2_reg_strength * ||coef||_2^2 for l2_reg_strength=1 and solve with
    # optimizer. Note that the problem is well conditioned such that we get accurate
    # results.
    l2_reg_strength = 1
    fun = partial(
        linear_loss.loss,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    grad = partial(
        linear_loss.gradient,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    coef_penalized_with_intercept = _special_minimize(
        fun, grad, coef_unpenalized, tol_NM=1e-6, tol=1e-14
    )

    linear_loss = LinearModelLoss(base_loss=model._get_loss(), fit_intercept=False)
    fun = partial(
        linear_loss.loss,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    grad = partial(
        linear_loss.gradient,
        X=X[:, :-1],
        y=y,
        sample_weight=sw,
        l2_reg_strength=l2_reg_strength,
    )
    coef_penalized_without_intercept = _special_minimize(
        fun, grad, coef_unpenalized[:-1], tol_NM=1e-6, tol=1e-14
    )

    # To be sure
    assert np.linalg.norm(coef_penalized_with_intercept) < np.linalg.norm(
        coef_unpenalized
    )

    return (
        model,
        X,
        y,
        coef_unpenalized,
        coef_penalized_with_intercept,
        coef_penalized_without_intercept,
        l2_reg_strength,
    )