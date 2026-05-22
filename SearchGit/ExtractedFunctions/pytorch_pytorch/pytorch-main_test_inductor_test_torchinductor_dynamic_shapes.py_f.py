        def f(x):
            y = x.item()
            return torch.ones(1, y, device=device) @ torch.ones(y, 1, device=device)