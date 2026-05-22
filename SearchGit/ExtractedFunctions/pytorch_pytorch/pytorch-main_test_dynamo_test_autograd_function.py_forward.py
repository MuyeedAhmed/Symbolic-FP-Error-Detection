            def forward(ctx, x, z):
                # Simple computation
                o = torch.matmul(x, x) @ x
                out = x.sin()
                # Mutate the nonlocal list
                z.append(out)
                return torch.cos(torch.sin(o)), torch.sin(x)