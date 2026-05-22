        def f(x, y):
            with torch.amp.autocast("cpu"):
                x = x @ y
            return x