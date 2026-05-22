    def fit(self, x_train: np.ndarray, y_train: np.ndarray) -> None:
        """
        Computes the polynomial regression model parameters using ordinary least squares
        (OLS) estimation:

        β = (XᵀX)⁻¹Xᵀy = X⁺y

        where X⁺ denotes the Moore-Penrose pseudoinverse of the design matrix X. This
        function computes X⁺ using singular value decomposition (SVD).

        References:
            - https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse
            - https://en.wikipedia.org/wiki/Singular_value_decomposition
            - https://en.wikipedia.org/wiki/Multicollinearity

        @param x_train: the predictor values x for model fitting
        @param y_train: the response values y for model fitting
        @raises ArithmeticError:    if X isn't full rank, then XᵀX is singular and β
                                    doesn't exist

        >>> x = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        >>> y = x**3 - 2 * x**2 + 3 * x - 5
        >>> poly_reg = PolynomialRegression(degree=3)
        >>> poly_reg.fit(x, y)
        >>> poly_reg.params
        array([-5.,  3., -2.,  1.])
        >>> poly_reg = PolynomialRegression(degree=20)
        >>> poly_reg.fit(x, y)
        Traceback (most recent call last):
        ...
        ArithmeticError: Design matrix is not full rank, can't compute coefficients

        Make sure errors don't grow too large:
        >>> coefs = np.array([-250, 50, -2, 36, 20, -12, 10, 2, -1, -15, 1])
        >>> y = PolynomialRegression._design_matrix(x, len(coefs) - 1) @ coefs
        >>> poly_reg = PolynomialRegression(degree=len(coefs) - 1)
        >>> poly_reg.fit(x, y)
        >>> np.allclose(poly_reg.params, coefs, atol=10e-3)
        True
        """
        X = PolynomialRegression._design_matrix(x_train, self.degree)  # noqa: N806
        _, cols = X.shape
        if np.linalg.matrix_rank(X) < cols:
            raise ArithmeticError(
                "Design matrix is not full rank, can't compute coefficients"
            )

        # np.linalg.pinv() computes the Moore-Penrose pseudoinverse using SVD
        self.params = np.linalg.pinv(X) @ y_train