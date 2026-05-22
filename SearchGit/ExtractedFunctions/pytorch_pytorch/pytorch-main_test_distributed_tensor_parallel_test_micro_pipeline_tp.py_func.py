        def func(inp: torch.Tensor, w1: torch.Tensor, w2: torch.Tensor) -> torch.Tensor:
            hidden = all_gather_tensor(inp, 0, (device_mesh, 0)) @ w1.t()
            full_hidden = all_gather_tensor(hidden, 0, (device_mesh, 1))
            full_hidden /= full_hidden.pow(2).sum().sqrt()
            hidden = reduce_scatter_tensor(full_hidden, "avg", 0, (device_mesh, 1))
            return reduce_scatter_tensor(hidden @ w2.t(), "avg", 0, (device_mesh, 0))