    def __init__(
        self,
        rng,
        n_samples=200,
        n_components=2,
        n_features=2,
        scale=50,
        dtype=np.float64,
    ):
        self.n_samples = n_samples
        self.n_components = n_components
        self.n_features = n_features

        self.weights = rng.rand(n_components).astype(dtype)
        self.weights = self.weights.astype(dtype) / self.weights.sum()
        self.means = rng.rand(n_components, n_features).astype(dtype) * scale
        self.covariances = {
            "spherical": 0.5 + rng.rand(n_components).astype(dtype),
            "diag": (0.5 + rng.rand(n_components, n_features).astype(dtype)) ** 2,
            "tied": make_spd_matrix(n_features, random_state=rng).astype(dtype),
            "full": np.array(
                [
                    make_spd_matrix(n_features, random_state=rng).astype(dtype) * 0.5
                    for _ in range(n_components)
                ]
            ),
        }
        self.precisions = {
            "spherical": 1.0 / self.covariances["spherical"],
            "diag": 1.0 / self.covariances["diag"],
            "tied": linalg.inv(self.covariances["tied"]),
            "full": np.array(
                [linalg.inv(covariance) for covariance in self.covariances["full"]]
            ),
        }

        self.X = dict(
            zip(
                COVARIANCE_TYPE,
                [
                    generate_data(
                        n_samples,
                        n_features,
                        self.weights,
                        self.means,
                        self.covariances,
                        covar_type,
                        dtype=dtype,
                    )
                    for covar_type in COVARIANCE_TYPE
                ],
            )
        )
        self.Y = np.hstack(
            [
                np.full(int(np.round(w * n_samples)), k, dtype=int)
                for k, w in enumerate(self.weights)
            ]
        )