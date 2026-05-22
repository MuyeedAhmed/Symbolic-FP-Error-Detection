        def mixed_inputs(x, y):
            out1 = x @ closure_weight1 + y
            out2 = x @ closure_weight2
            return (out1, out2)