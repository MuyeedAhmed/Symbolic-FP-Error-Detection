    def get_expected_eq_scales(self, model, x):
        """ For each module in the graph, we want to calculate the equalization
        scale at that point. This only works for models containing single or
        connected linear layers.
        """
        exp_eq_scales = []
        for _, module in model.named_children():
            weight = module.weight.detach().numpy()
            bias = module.bias.detach().numpy()

            eq_scale = self.calculate_equalization_scale_ref(x, weight)
            exp_eq_scales.append(eq_scale)

            x = x @ weight.T + bias

        return exp_eq_scales