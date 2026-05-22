        def f(x, w):
            with torch.inference_mode():
                return (x @ w).sum()