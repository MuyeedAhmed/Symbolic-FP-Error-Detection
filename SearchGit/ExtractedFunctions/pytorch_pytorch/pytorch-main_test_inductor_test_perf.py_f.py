        def f(a, b):
            a = a @ a
            return torch.constant_pad_nd(torch.cat([a, b]), [2, 2], 0.5)