def _omp_path_residues(
    X_train,
    y_train,
    X_test,
    y_test,
    copy=True,
    fit_intercept=True,
    max_iter=100,
):
    """Compute the residues on left-out data for a full LARS path.

    Parameters
    ----------
    X_train : ndarray of shape (n_samples, n_features)
        The data to fit the LARS on.

    y_train : ndarray of shape (n_samples)
        The target variable to fit LARS on.

    X_test : ndarray of shape (n_samples, n_features)
        The data to compute the residues on.

    y_test : ndarray of shape (n_samples)
        The target variable to compute the residues on.

    copy : bool, default=True
        Whether X_train, X_test, y_train and y_test should be copied.  If
        False, they may be overwritten.

    fit_intercept : bool, default=True
        Whether to calculate the intercept for this model. If set
        to false, no intercept will be used in calculations
        (i.e. data is expected to be centered).

    max_iter : int, default=100
        Maximum numbers of iterations to perform, therefore maximum features
        to include. 100 by default.

    Returns
    -------
    residues : ndarray of shape (n_samples, max_features)
        Residues of the prediction on the test data.
    """

    if copy:
        X_train = X_train.copy()
        y_train = y_train.copy()
        X_test = X_test.copy()
        y_test = y_test.copy()

    if fit_intercept:
        X_mean = X_train.mean(axis=0)
        X_train -= X_mean
        X_test -= X_mean
        y_mean = y_train.mean(axis=0)
        y_train = as_float_array(y_train, copy=False)
        y_train -= y_mean
        y_test = as_float_array(y_test, copy=False)
        y_test -= y_mean

    coefs = orthogonal_mp(
        X_train,
        y_train,
        n_nonzero_coefs=max_iter,
        tol=None,
        precompute=False,
        copy_X=False,
        return_path=True,
    )
    if coefs.ndim == 1:
        coefs = coefs[:, np.newaxis]

    return np.dot(coefs.T, X_test.T) - y_test