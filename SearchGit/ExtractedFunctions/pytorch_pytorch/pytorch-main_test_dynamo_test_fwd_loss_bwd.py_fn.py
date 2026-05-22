        def fn(a, b):
            loss = (a @ b).sum()
            loss.backward()
            return a.grad