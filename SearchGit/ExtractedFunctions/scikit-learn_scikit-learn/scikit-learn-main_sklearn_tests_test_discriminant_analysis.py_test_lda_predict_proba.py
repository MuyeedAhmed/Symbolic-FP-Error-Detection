def test_lda_predict_proba(solver, n_classes):
    def generate_dataset(n_samples, centers, covariances, random_state=None):
        """Generate a multivariate normal data given some centers and
        covariances"""
        rng = check_random_state(random_state)
        X = np.vstack(
            [
                rng.multivariate_normal(mean, cov, size=n_samples // len(centers))
                for mean, cov in zip(centers, covariances)
            ]
        )
        y = np.hstack(
            [[clazz] * (n_samples // len(centers)) for clazz in range(len(centers))]
        )
        return X, y

    blob_centers = np.array([[0, 0], [-10, 40], [-30, 30]])[:n_classes]
    blob_stds = np.array([[[10, 10], [10, 100]]] * len(blob_centers))
    X, y = generate_dataset(
        n_samples=90000, centers=blob_centers, covariances=blob_stds, random_state=42
    )
    lda = LinearDiscriminantAnalysis(
        solver=solver, store_covariance=True, shrinkage=None
    ).fit(X, y)
    # check that the empirical means and covariances are close enough to the
    # one used to generate the data
    assert_allclose(lda.means_, blob_centers, atol=1e-1)
    assert_allclose(lda.covariance_, blob_stds[0], atol=1)

    # implement the method to compute the probability given in The Elements
    # of Statistical Learning (cf. p.127, Sect. 4.4.5 "Logistic Regression
    # or LDA?")
    precision = linalg.inv(blob_stds[0])
    alpha_k = []
    alpha_k_0 = []
    for clazz in range(len(blob_centers) - 1):
        alpha_k.append(
            np.dot(precision, (blob_centers[clazz] - blob_centers[-1])[:, np.newaxis])
        )
        alpha_k_0.append(
            np.dot(
                -0.5 * (blob_centers[clazz] + blob_centers[-1])[np.newaxis, :],
                alpha_k[-1],
            )
        )

    sample = np.array([[-22, 22]])

    def discriminant_func(sample, coef, intercept, clazz):
        return np.exp(intercept[clazz] + np.dot(sample, coef[clazz])).item()

    prob = np.array(
        [
            float(
                discriminant_func(sample, alpha_k, alpha_k_0, clazz)
                / (
                    1
                    + sum(
                        [
                            discriminant_func(sample, alpha_k, alpha_k_0, clazz)
                            for clazz in range(n_classes - 1)
                        ]
                    )
                )
            )
            for clazz in range(n_classes - 1)
        ]
    )

    prob_ref = 1 - np.sum(prob)

    # check the consistency of the computed probability
    # all probabilities should sum to one
    prob_ref_2 = float(
        1
        / (
            1
            + sum(
                [
                    discriminant_func(sample, alpha_k, alpha_k_0, clazz)
                    for clazz in range(n_classes - 1)
                ]
            )
        )
    )

    assert prob_ref == pytest.approx(prob_ref_2)
    # check that the probability of LDA are close to the theoretical
    # probabilities
    assert_allclose(
        lda.predict_proba(sample), np.hstack([prob, prob_ref])[np.newaxis], atol=1e-2
    )