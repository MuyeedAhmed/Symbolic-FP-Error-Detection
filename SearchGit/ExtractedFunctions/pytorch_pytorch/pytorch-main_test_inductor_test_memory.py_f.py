        def f(a, b, c):
            return (a @ c).sum(dim=-1) + (b @ c).sum(dim=-1)