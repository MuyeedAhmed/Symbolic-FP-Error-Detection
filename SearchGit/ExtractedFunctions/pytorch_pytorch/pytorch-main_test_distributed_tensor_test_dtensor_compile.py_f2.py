        def f2(x, y):
            torch._check(x.size(0) <= 16384)
            torch._check(y.size(1) <= 32)
            return x @ y