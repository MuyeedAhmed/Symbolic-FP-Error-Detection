        def fn(x, w):
            y = x @ w
            y = y.view(2, 4, 16).permute(1, 0, 2).permute(1, 0, 2).reshape(8, 16)
            return y + x