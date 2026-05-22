    def predict(self, X):
        """Perform classification on an array of test vectors X.

        Parameters
        ----------
        X : array-like of shape (n_samples, n_features) or list of object
            Query points where the GP is evaluated for classification.

        Returns
        -------
        C : ndarray of shape (n_samples,)
            Predicted target values for X, values are from ``classes_``
        """
        check_is_fitted(self)

        # As discussed on Section 3.4.2 of GPML, for making hard binary
        # decisions, it is enough to compute the MAP of the posterior and
        # pass it through the link function
        K_star = self.kernel_(self.X_train_, X)  # K_star =k(x_star)
        f_star = K_star.T.dot(self.y_train_ - self.pi_)  # Algorithm 3.2,Line 4

        return np.where(f_star > 0, self.classes_[1], self.classes_[0])