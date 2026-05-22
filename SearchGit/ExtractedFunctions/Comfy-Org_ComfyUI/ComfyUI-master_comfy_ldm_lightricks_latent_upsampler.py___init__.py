    def __init__(self, dims: int, stride: int):
        super().__init__()
        assert dims in (2, 3)
        assert stride >= 1 and isinstance(stride, int)
        self.dims = dims
        self.stride = stride

        # 5x5 separable binomial kernel [1,4,6,4,1] (outer product), normalized
        k = torch.tensor([1.0, 4.0, 6.0, 4.0, 1.0])
        k2d = k[:, None] @ k[None, :]
        k2d = (k2d / k2d.sum()).float()  # shape (5,5)
        self.register_buffer("kernel", k2d[None, None, :, :])  # (1,1,5,5)