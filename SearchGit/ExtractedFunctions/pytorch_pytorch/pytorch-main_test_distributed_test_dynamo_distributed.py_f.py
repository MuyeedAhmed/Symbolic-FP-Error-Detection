            def f(a, b, c):
                res = (
                    torch.sum((a @ b) + 1.0)
                    + torch.sum(torch.relu(b @ c))
                    + torch.sum(c @ a)
                )

                return res