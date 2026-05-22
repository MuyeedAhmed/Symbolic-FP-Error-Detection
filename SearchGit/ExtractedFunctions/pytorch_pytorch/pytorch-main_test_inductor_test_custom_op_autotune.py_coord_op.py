        def coord_op(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            if x.shape[0] == 128:
                return torch.empty(
                    x.shape[0],
                    weight.shape[1],
                    dtype=weight.dtype,
                    device=weight.device,
                )
            return x @ weight