        def f(idx, x):
            x = x.select(0, idx.item())
            return x @ x