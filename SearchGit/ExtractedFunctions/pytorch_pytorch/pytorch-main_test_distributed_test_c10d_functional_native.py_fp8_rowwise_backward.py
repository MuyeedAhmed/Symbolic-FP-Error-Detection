        def fp8_rowwise_backward(in_, w, out_grad):
            out_grad_fp8, scale_out_grad = scale(out_grad)
            w_fp8, scale_w = scale(w.t().contiguous())
            out_grad_fp8 = funcol.all_gather_tensor(
                out_grad_fp8, gather_dim=0, group=torch.distributed.group.WORLD
            )
            scale_out_grad = funcol.all_gather_tensor(
                scale_out_grad, gather_dim=0, group=torch.distributed.group.WORLD
            )
            in_grad = torch._scaled_mm(
                out_grad_fp8,
                w_fp8.t(),
                scale_a=scale_out_grad,
                scale_b=scale_w.t(),
                out_dtype=torch.bfloat16,
            )

            out_grad = funcol.all_gather_tensor(
                out_grad.t().contiguous(),
                gather_dim=0,
                group=torch.distributed.group.WORLD,
            )
            w_grad = out_grad @ in_

            return in_grad, w_grad