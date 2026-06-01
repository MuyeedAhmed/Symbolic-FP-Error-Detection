        def fn(args):
            (x,) = args
            return torch.linalg.inv(x)