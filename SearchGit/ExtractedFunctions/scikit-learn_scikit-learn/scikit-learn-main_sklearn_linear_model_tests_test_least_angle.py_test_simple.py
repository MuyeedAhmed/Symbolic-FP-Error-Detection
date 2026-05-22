def test_simple():
    # Principle of Lars is to keep covariances tied and decreasing

    # also test verbose output
    import sys
    from io import StringIO

    old_stdout = sys.stdout
    try:
        sys.stdout = StringIO()

        _, _, coef_path_ = linear_model.lars_path(X, y, method="lar", verbose=10)

        sys.stdout = old_stdout

        for i, coef_ in enumerate(coef_path_.T):
            res = y - np.dot(X, coef_)
            cov = np.dot(X.T, res)
            C = np.max(abs(cov))
            eps = 1e-3
            ocur = len(cov[C - eps < abs(cov)])
            if i < X.shape[1]:
                assert ocur == i + 1
            else:
                # no more than max_pred variables can go into the active set
                assert ocur == X.shape[1]
    finally:
        sys.stdout = old_stdout