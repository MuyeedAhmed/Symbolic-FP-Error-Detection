    def back_propagation(self) -> None:
        """
        Function for fine-tuning the weights of the neural net based on the
        error rate obtained in the previous epoch (i.e., iteration).
        Updation is done using derivative of sogmoid activation function.

        >>> input_val = np.array(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
        >>> output_val = np.array(([0], [0], [0]), dtype=float)
        >>> nn = TwoHiddenLayerNeuralNetwork(input_val, output_val)
        >>> res = nn.feedforward()
        >>> nn.back_propagation()
        >>> updated_weights = nn.second_hidden_layer_and_output_layer_weights
        >>> bool((res == updated_weights).all())
        False
        """

        updated_second_hidden_layer_and_output_layer_weights = np.dot(
            self.layer_between_first_hidden_layer_and_second_hidden_layer.T,
            2
            * (self.output_array - self.predicted_output)
            * sigmoid_derivative(self.predicted_output),
        )
        updated_first_hidden_layer_and_second_hidden_layer_weights = np.dot(
            self.layer_between_input_and_first_hidden_layer.T,
            np.dot(
                2
                * (self.output_array - self.predicted_output)
                * sigmoid_derivative(self.predicted_output),
                self.second_hidden_layer_and_output_layer_weights.T,
            )
            * sigmoid_derivative(
                self.layer_between_first_hidden_layer_and_second_hidden_layer
            ),
        )
        updated_input_layer_and_first_hidden_layer_weights = np.dot(
            self.input_array.T,
            np.dot(
                np.dot(
                    2
                    * (self.output_array - self.predicted_output)
                    * sigmoid_derivative(self.predicted_output),
                    self.second_hidden_layer_and_output_layer_weights.T,
                )
                * sigmoid_derivative(
                    self.layer_between_first_hidden_layer_and_second_hidden_layer
                ),
                self.first_hidden_layer_and_second_hidden_layer_weights.T,
            )
            * sigmoid_derivative(self.layer_between_input_and_first_hidden_layer),
        )

        self.input_layer_and_first_hidden_layer_weights += (
            updated_input_layer_and_first_hidden_layer_weights
        )
        self.first_hidden_layer_and_second_hidden_layer_weights += (
            updated_first_hidden_layer_and_second_hidden_layer_weights
        )
        self.second_hidden_layer_and_output_layer_weights += (
            updated_second_hidden_layer_and_output_layer_weights
        )