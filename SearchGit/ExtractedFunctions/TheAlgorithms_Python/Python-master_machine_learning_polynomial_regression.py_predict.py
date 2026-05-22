    def predict(self, data: np.ndarray) -> np.ndarray:
        """
        Computes the predicted response values y for the given input data by
        constructing the design matrix X and evaluating y = Xβ.

        @param data:    the predictor values x for prediction
        @returns:       the predicted response values y = Xβ
        @raises ArithmeticError:    if this function is called before the model
                                    parameters are fit

        >>> x = np.array([0, 1, 2, 3, 4])
        >>> y = x**3 - 2 * x**2 + 3 * x - 5
        >>> poly_reg = PolynomialRegression(degree=3)
        >>> poly_reg.fit(x, y)
        >>> poly_reg.predict(np.array([-1]))
        array([-11.])
        >>> poly_reg.predict(np.array([-2]))
        array([-27.])
        >>> poly_reg.predict(np.array([6]))
        array([157.])
        >>> PolynomialRegression(degree=3).predict(x)
        Traceback (most recent call last):
        ...
        ArithmeticError: Predictor hasn't been fit yet
        """
        if self.params is None:
            raise ArithmeticError("Predictor hasn't been fit yet")

        return PolynomialRegression._design_matrix(data, self.degree) @ self.params