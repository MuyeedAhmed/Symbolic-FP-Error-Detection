        def impl_mm_bad_schema(x: torch.Tensor, y: torch.Tensor) -> torch.Tensor:
            return x @ y - 1