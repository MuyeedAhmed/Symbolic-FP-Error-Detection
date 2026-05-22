        def fn(x):
            z = torch.ones(3, device=device, dtype=dtype)
            return grad(lambda x: z @ x)(x)