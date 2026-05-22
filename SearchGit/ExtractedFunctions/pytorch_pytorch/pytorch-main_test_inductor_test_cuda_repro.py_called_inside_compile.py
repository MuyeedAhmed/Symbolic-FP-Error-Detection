        def called_inside_compile(x, w, b):
            a = x @ w + b
            return torch.sigmoid(a)