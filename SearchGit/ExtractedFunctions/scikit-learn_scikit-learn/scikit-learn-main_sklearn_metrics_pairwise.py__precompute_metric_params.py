def _precompute_metric_params(X, Y, metric=None, **kwds):
    """Precompute data-derived metric parameters if not provided."""
    if metric == "seuclidean" and "V" not in kwds:
        if X is Y:
            V = np.var(X, axis=0, ddof=1)
        else:
            raise ValueError(
                "The 'V' parameter is required for the seuclidean metric "
                "when Y is passed."
            )
        return {"V": V}
    if metric == "mahalanobis" and "VI" not in kwds:
        if X is Y:
            VI = np.linalg.inv(np.cov(X.T)).T
        else:
            raise ValueError(
                "The 'VI' parameter is required for the mahalanobis metric "
                "when Y is passed."
            )
        return {"VI": VI}
    return {}