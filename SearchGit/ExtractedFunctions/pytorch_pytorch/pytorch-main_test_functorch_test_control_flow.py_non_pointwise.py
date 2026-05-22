    def non_pointwise(x: torch.Tensor, y: torch.Tensor):
        W = torch.arange(4, dtype=torch.float, device=x.device).view(2, 2)
        return x @ W + y @ W