        def fn(x):
            return torch.linalg.inv(torch.linalg.cholesky(x))