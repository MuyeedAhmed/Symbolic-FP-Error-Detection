        def f(x):
            return (x @ x.transpose(0, 1)).relu()