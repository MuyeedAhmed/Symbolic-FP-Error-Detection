            def forward(ctx, x, w):
                # x: [M, N]
                # w: [N, K]
                ctx.save_for_backward(x, w)
                return x @ w