        def f(a, b):
            return (a @ b).to(torch.float32).sum(dim=1)