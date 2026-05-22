        def cudagraph_patch_op(x: torch.Tensor, weight: torch.Tensor) -> torch.Tensor:
            return x @ weight