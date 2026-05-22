    def h(self, x: torch.Tensor, base_out: torch.Tensor) -> torch.Tensor:
        """
        Additive bypass component for LoHa: h(x) = diff_weight @ x

        WARNING: Inefficient - computes full Hadamard product each forward.

        Note:
            Does not access original model weights - bypass mode is designed
            for quantized models where weights may not be accessible.

        Args:
            x: Input tensor
            base_out: Output from base forward (unused, for API consistency)

        Reference: LyCORIS functional/loha.py bypass_forward_diff
        """
        _warn_loha_bypass_inefficient()

        # FUNC_LIST: [None, None, F.linear, F.conv1d, F.conv2d, F.conv3d]
        FUNC_LIST = [None, None, F.linear, F.conv1d, F.conv2d, F.conv3d]

        v = self.weights
        # v[0]=w1a, v[1]=w1b, v[2]=alpha, v[3]=w2a, v[4]=w2b, v[5]=t1, v[6]=t2, v[7]=dora
        w1a = v[0]
        w1b = v[1]
        alpha = v[2]
        w2a = v[3]
        w2b = v[4]
        t1 = v[5]
        t2 = v[6]

        # Compute scale
        rank = w1b.shape[0]
        scale = (alpha / rank if alpha is not None else 1.0) * getattr(
            self, "multiplier", 1.0
        )

        # Cast dtype
        w1a = w1a.to(dtype=x.dtype)
        w1b = w1b.to(dtype=x.dtype)
        w2a = w2a.to(dtype=x.dtype)
        w2b = w2b.to(dtype=x.dtype)

        # Use module info from bypass injection, not weight dimension
        is_conv = getattr(self, "is_conv", False)
        conv_dim = getattr(self, "conv_dim", 0)
        kw_dict = getattr(self, "kw_dict", {})

        # Compute diff weight using Hadamard product
        if t1 is not None and t2 is not None:
            t1 = t1.to(dtype=x.dtype)
            t2 = t2.to(dtype=x.dtype)
            m1 = torch.einsum("i j k l, j r, i p -> p r k l", t1, w1b, w1a)
            m2 = torch.einsum("i j k l, j r, i p -> p r k l", t2, w2b, w2a)
            diff_weight = (m1 * m2) * scale
        else:
            m1 = w1a @ w1b
            m2 = w2a @ w2b
            diff_weight = (m1 * m2) * scale

        if is_conv:
            op = FUNC_LIST[conv_dim + 2]
            kernel_size = getattr(self, "kernel_size", (1,) * conv_dim)
            in_channels = getattr(self, "in_channels", None)

            # Reshape 2D diff_weight to conv format using kernel_size
            # diff_weight: [out_channels, in_channels * prod(kernel_size)] -> [out_channels, in_channels, *kernel_size]
            if diff_weight.dim() == 2:
                if in_channels is not None:
                    diff_weight = diff_weight.view(
                        diff_weight.shape[0], in_channels, *kernel_size
                    )
                else:
                    diff_weight = diff_weight.view(
                        *diff_weight.shape, *([1] * conv_dim)
                    )
        else:
            op = F.linear
            kw_dict = {}

        return op(x, diff_weight, **kw_dict)