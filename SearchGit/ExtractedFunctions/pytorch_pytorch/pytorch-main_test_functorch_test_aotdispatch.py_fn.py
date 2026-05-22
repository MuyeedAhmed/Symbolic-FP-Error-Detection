        def fn(p, x):
            y = x @ x
            y.add_(2)
            return (x.sum() + y.view(1, 4).sum(),)