    def backward(ctx, dY: torch.Tensor):
        W, W_quant, S = ctx.custom_saved_tensors
        A, B, X = ctx.saved_tensors

        batch, seq_len, hd = X.shape
        dY = dY.reshape(-1, dY.shape[-1])  # Must be reshape
        X = X.reshape(-1, X.shape[-1])  # Must be reshape
        dtype = X.dtype

        A, B = A.to(dtype), B.to(dtype)

        A, B = A.t(), B.t()

        d_A = torch.empty_like(A)
        d_B = torch.empty_like(B)

        ### Weight projection LoRA weights
        # Weight projection
        # d_A = X.t() @ (dY @ B.t())
        # d_B = (A.t() @ X.t()) @ dY
        # d_A *= S
        # d_B *= S
        d_A.addmm_(X.t(), dY @ B.t(), alpha = S, beta = 0)
        d_B.addmm_(A.t() @ X.t(), dY, alpha = S, beta = 0)

        # Get derivative for dX
        W = fast_dequantize(W.t(), W_quant)
        dX = dY @ W.t()
        del W
        # dX += dY @ B.to(dtype).t() @ (S * A.to(dtype).t())
        dX.addmm_(dY @ B.t(), A.t(), alpha = S)

        # W, W_quant, A, B, S
        return dX.view(batch, seq_len, hd), None, None, d_A.t(), d_B.t(), None