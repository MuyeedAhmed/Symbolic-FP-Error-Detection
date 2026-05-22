            def backward(ctx, g_out):
                # g_out: [M, K]
                x, w = ctx.saved_tensors
                g_x = g_out @ w.T
                g_w = x.T @ g_out
                w.main_grad = g_w.float()
                return g_x, None