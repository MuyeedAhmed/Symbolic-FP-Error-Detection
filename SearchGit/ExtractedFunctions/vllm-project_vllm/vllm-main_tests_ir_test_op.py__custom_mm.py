    def _custom_mm(
        x: torch.Tensor, y: torch.Tensor, bias: torch.Tensor | None = None
    ) -> torch.Tensor:
        tmp = x @ y
        return tmp if bias is None else tmp + bias