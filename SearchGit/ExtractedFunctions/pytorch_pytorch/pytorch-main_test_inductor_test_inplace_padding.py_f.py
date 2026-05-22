        def f(x, y):
            x = aten.constant_pad_nd(x, (0, 8, 0, 0), 12345.0)
            x.add_(1)
            return x @ y