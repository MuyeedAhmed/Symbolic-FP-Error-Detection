    def get_expected_weights_bias(self, model, x, exp_eq_scales):
        """ For each module in the graph, we want to calculate the expected
        scaled weight and bias values. This only works for models containing
        single or connected linear layers.
        """
        exp_weights = []
        exp_bias = []
        for i, (_, module) in enumerate(model.named_children()):
            weight = module.weight.detach().numpy()
            bias = module.bias.detach().numpy()

            scaled_weight = weight * np.reciprocal(exp_eq_scales[i])
            scaled_bias = bias
            if i + 1 < len(exp_eq_scales):
                scaled_weight = (scaled_weight.T * exp_eq_scales[i + 1]).T
                scaled_bias = (scaled_bias.T * exp_eq_scales[i + 1]).T

            exp_weights.append(scaled_weight)
            exp_bias.append(scaled_bias)

            x = x @ weight.T + bias

        return exp_weights, exp_bias