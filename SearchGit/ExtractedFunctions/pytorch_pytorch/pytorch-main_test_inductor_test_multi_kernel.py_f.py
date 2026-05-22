        def f(x, y):
            return x.sum(dim=-1, keepdims=True) * (y @ y)