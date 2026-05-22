        def grad_in_loop(x, y):
            for _ in range(100):
                x = y @ x
            return x