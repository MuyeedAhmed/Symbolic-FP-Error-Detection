        def fn(w_0, w_1, w_2, x):
            x = (x @ w_0).relu()
            x = (x @ w_1).relu()
            x = (x @ w_2).relu()
            return x