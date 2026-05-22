        def gn(x, z):
            o = torch.matmul(x, x) @ x
            out = x.sin()
            z.append(out)
            return torch.cos(torch.sin(o)), torch.sin(x)