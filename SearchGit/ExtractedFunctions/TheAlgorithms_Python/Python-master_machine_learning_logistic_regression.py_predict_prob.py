    def predict_prob(x):
        return sigmoid_function(
            np.dot(x, theta)
        )  # predicting the value of probability from the logistic regression algorithm