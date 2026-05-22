    def predict(self, input_arr: np.ndarray) -> int:
        """
        Predict's the output for the given input values using
        the trained neural network.

        The output value given by the model ranges in-between 0 and 1.
        The predict function returns 1 if the model value is greater
        than the threshold value else returns 0,
        as the real output values are in binary.

        >>> input_val = np.array(([0, 0, 0], [0, 1, 0], [0, 0, 1]), dtype=float)
        >>> output_val = np.array(([0], [1], [1]), dtype=float)
        >>> nn = TwoHiddenLayerNeuralNetwork(input_val, output_val)
        >>> nn.train(output_val, 1000, False)
        >>> nn.predict([0, 1, 0]) in (0, 1)
        True
        """

        # Input values for which the predictions are to be made.
        self.array = input_arr

        self.layer_between_input_and_first_hidden_layer = sigmoid(
            np.dot(self.array, self.input_layer_and_first_hidden_layer_weights)
        )

        self.layer_between_first_hidden_layer_and_second_hidden_layer = sigmoid(
            np.dot(
                self.layer_between_input_and_first_hidden_layer,
                self.first_hidden_layer_and_second_hidden_layer_weights,
            )
        )

        self.layer_between_second_hidden_layer_and_output = sigmoid(
            np.dot(
                self.layer_between_first_hidden_layer_and_second_hidden_layer,
                self.second_hidden_layer_and_output_layer_weights,
            )
        )

        return int((self.layer_between_second_hidden_layer_and_output > 0.6)[0])