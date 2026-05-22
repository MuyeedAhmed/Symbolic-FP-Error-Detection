    def feedforward(self) -> np.ndarray:
        """
        The information moves in only one direction i.e. forward from the input nodes,
        through the two hidden nodes and to the output nodes.
        There are no cycles or loops in the network.

        Return layer_between_second_hidden_layer_and_output
            (i.e the last layer of the neural network).

        >>> input_val = np.array(([0, 0, 0], [0, 0, 0], [0, 0, 0]), dtype=float)
        >>> output_val = np.array(([0], [0], [0]), dtype=float)
        >>> nn = TwoHiddenLayerNeuralNetwork(input_val, output_val)
        >>> res = nn.feedforward()
        >>> array_sum = np.sum(res)
        >>> bool(np.isnan(array_sum))
        False
        """
        # Layer_between_input_and_first_hidden_layer is the layer connecting the
        # input nodes with the first hidden layer nodes.
        self.layer_between_input_and_first_hidden_layer = sigmoid(
            np.dot(self.input_array, self.input_layer_and_first_hidden_layer_weights)
        )

        # layer_between_first_hidden_layer_and_second_hidden_layer is the layer
        # connecting the first hidden set of nodes with the second hidden set of nodes.
        self.layer_between_first_hidden_layer_and_second_hidden_layer = sigmoid(
            np.dot(
                self.layer_between_input_and_first_hidden_layer,
                self.first_hidden_layer_and_second_hidden_layer_weights,
            )
        )

        # layer_between_second_hidden_layer_and_output is the layer connecting
        # second hidden layer with the output node.
        self.layer_between_second_hidden_layer_and_output = sigmoid(
            np.dot(
                self.layer_between_first_hidden_layer_and_second_hidden_layer,
                self.second_hidden_layer_and_output_layer_weights,
            )
        )

        return self.layer_between_second_hidden_layer_and_output