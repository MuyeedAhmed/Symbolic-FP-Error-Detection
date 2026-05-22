        def loss(A, x1, x2):
            x2_hat = (A @ (x1.T)).T
            res = x2 - x2_hat
            res_sqr = res**2
            return res_sqr.sum()