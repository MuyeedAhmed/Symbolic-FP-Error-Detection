        def compute_loss(weight, x, t):
            y = x @ weight
            return ((y - t) ** 2).sum()