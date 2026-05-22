        def add(carry, x):
            ret = (carry @ param_buffer[0]["param"]) @ x + param_buffer[1][0]
            return ret, ret.sum()