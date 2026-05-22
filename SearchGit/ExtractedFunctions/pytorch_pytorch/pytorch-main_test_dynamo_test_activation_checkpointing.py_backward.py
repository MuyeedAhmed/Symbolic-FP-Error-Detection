            def backward(ctx, grad):
                x, w = ctx.saved_tensors
                return grad @ w.t(), x.t() @ grad