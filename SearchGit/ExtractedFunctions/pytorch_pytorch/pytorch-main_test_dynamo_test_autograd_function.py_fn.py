        def fn(w):
            y = torch.ones(2, 2) @ w
            TimesThreeInplace.apply(y)
            return y.sum()