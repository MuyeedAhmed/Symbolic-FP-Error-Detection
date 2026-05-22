def test_rank_deficient_design():
    # consistency test that checks that LARS Lasso is handling rank
    # deficient input data (with n_features < rank) in the same way
    # as coordinate descent Lasso
    y = [5, 0, 5]
    for X in ([[5, 0], [0, 5], [10, 10]], [[10, 10, 0], [1e-32, 0, 0], [0, 0, 1]]):
        # To be able to use the coefs to compute the objective function,
        # we need to turn off normalization
        lars = linear_model.LassoLars(0.1)
        coef_lars_ = lars.fit(X, y).coef_
        obj_lars = 1.0 / (2.0 * 3.0) * linalg.norm(
            y - np.dot(X, coef_lars_)
        ) ** 2 + 0.1 * linalg.norm(coef_lars_, 1)
        coord_descent = linear_model.Lasso(0.1, tol=1e-6)
        coef_cd_ = coord_descent.fit(X, y).coef_
        obj_cd = (1.0 / (2.0 * 3.0)) * linalg.norm(
            y - np.dot(X, coef_cd_)
        ) ** 2 + 0.1 * linalg.norm(coef_cd_, 1)
        assert obj_lars < obj_cd * (1.0 + 1e-8)