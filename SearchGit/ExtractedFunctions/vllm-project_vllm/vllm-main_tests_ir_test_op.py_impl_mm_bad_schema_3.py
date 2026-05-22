        def impl_mm_bad_schema_3(
            x: torch.Tensor, y: torch.Tensor, bias: torch.Tensor
        ) -> torch.Tensor:
            return x @ y + bias - 5