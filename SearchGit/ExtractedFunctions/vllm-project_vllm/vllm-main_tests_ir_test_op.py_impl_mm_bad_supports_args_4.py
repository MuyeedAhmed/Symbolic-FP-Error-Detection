        def impl_mm_bad_supports_args_4(
            x: torch.Tensor, y: torch.Tensor, bias: torch.Tensor | None = None
        ) -> torch.Tensor:
            return x @ y + 30