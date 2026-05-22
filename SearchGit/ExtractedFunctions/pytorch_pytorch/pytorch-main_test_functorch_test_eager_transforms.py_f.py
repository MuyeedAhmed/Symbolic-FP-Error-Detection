        def f(x, y, z):
            return ((x.relu() * x) @ y.sin() @ z).sum()