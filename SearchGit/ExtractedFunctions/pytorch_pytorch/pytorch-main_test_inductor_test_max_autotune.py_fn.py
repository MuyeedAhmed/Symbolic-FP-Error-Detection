        def fn(a, b, c):
            a = (a @ b) @ c
            a, b, c = (t.to(torch.float16) for t in [a, b, c])
            return (a @ b) @ c