    def get_expected_inp_act_vals(self, model, x, exp_eq_scales, exp_weights, exp_bias):
        """ For each module in the graph, we want to calculate the expected
        min/max values for every input activation node. This only works for
        models containing only single or connected linear layers.
        """
        x = x * exp_eq_scales[0]

        exp_inp_activation_vals = []
        for i, _ in enumerate(model.named_children()):
            exp_inp_activation_vals.append((x.min(), x.max()))
            x = x @ exp_weights[i].T + exp_bias[i]

        exp_inp_activation_vals.append((x.min(), x.max()))
        return exp_inp_activation_vals