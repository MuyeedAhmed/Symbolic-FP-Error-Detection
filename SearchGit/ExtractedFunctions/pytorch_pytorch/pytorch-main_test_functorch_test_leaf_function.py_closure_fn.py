        def closure_fn(x, y):
            if config["use_double"]:
                return (x @ y * 2,)
            return (x @ y * 3,)