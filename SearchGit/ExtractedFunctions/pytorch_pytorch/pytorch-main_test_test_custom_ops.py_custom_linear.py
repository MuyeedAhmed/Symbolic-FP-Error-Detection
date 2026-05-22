        def custom_linear(x: Tensor, weight: Tensor, bias: Tensor) -> Tensor:
            nonlocal called_impl
            called_impl = True
            x_np = x.numpy()
            w_np = weight.numpy()
            b_np = bias.numpy()
            out_np = np.add(x_np @ w_np.T, bias)
            return out_np