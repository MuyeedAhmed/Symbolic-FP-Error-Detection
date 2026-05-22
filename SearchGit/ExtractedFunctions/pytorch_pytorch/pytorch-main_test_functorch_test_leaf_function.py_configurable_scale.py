        def configurable_scale(x, y):
            if closure_config["use_double_scale"]:
                return (x @ y * 2,)
            else:
                return (x @ y * 3,)