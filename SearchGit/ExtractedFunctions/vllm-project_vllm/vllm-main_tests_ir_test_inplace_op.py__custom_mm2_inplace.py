def _custom_mm2_inplace(x: Tensor, w: Tensor) -> Tensor:
    x.copy_(x @ w + 2)
    return x