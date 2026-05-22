        def impl_mm_bad_schema_2(
            x: torch.Tensor, y: torch.Tensor, b: torch.Tensor | None = None
        ) -> torch.Tensor:
            return x @ y + b - 2