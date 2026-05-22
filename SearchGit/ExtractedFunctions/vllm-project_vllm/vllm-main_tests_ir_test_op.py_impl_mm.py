    def impl_mm(
        x: torch.Tensor, y: torch.Tensor, bias: torch.Tensor | None = None
    ) -> torch.Tensor:
        tmp = x @ y
        return tmp + 50 if bias is None else tmp + bias + 100