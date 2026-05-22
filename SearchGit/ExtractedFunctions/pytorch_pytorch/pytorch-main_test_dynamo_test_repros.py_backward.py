            def backward(ctx, grad_out):
                a, b_t, amax_b_t = ctx.saved_tensors
                # Workaround for https://github.com/pytorch/pytorch/issues/141881.
                # The partitioner would pre-compute the transposed scaling of the weight
                # in the forward (as it's most efficient, but it actually uses too much
                # memory). We prevent that by making the scaling depend on the gradient
                # in a way that has no effect and will be optimized away later.
                # Care is needed to support tensor parallelism and circumvent bugs.
                #        b_t = b_t + grad_out[:1, :, None].squeeze(0) * 0
                if ctx.a_requires_grad:
                    b = b_t.t().contiguous()
                    amax_grad_out = grad_out.abs().unflatten(-1, (1, -1)).amax(dim=-1)
                    amax_b = amax_b_t.t().unflatten(-1, (1, -1)).amax(dim=-1)
                    amax_b = amax_b.repeat_interleave(
                        b.shape[0] // amax_b.shape[0], dim=0, output_size=b.shape[0]
                    )
                    grad_a = matmul(grad_out, amax_grad_out, b, amax_b, None)
                else:
                    grad_a = None
                if ctx.b_requires_grad:
                    grad_b = grad_out.t() @ a
                else:
                    grad_b = None
                if ctx.bias_requires_grad:
                    grad_bias = grad_out.sum(dim=0)
                else:
                    grad_bias = None
                return grad_a, grad_b, grad_bias