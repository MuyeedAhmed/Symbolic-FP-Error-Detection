        def fn2(x, y):
            yy = y @ y
            return x * 2 + yy.view(25)