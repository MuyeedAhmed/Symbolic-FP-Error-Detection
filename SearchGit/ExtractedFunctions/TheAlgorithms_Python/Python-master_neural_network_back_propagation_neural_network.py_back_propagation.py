    def back_propagation(self, gradient):
        gradient_activation = self.cal_gradient()  # i * i 维
        gradient = np.asmatrix(np.dot(gradient.T, gradient_activation))

        self._gradient_weight = np.asmatrix(self.xdata)
        self._gradient_bias = -1
        self._gradient_x = self.weight

        self.gradient_weight = np.dot(gradient.T, self._gradient_weight.T)
        self.gradient_bias = gradient * self._gradient_bias
        self.gradient = np.dot(gradient, self._gradient_x).T
        # upgrade: the Negative gradient direction
        self.weight = self.weight - self.learn_rate * self.gradient_weight
        self.bias = self.bias - self.learn_rate * self.gradient_bias.T
        # updates the weights and bias according to learning rate (0.3 if undefined)
        return self.gradient