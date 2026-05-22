        def f1(x, y):
            torch._check(x.size(0) <= 32)
            torch._check(y.size(1) <= 16384)
            return x @ y