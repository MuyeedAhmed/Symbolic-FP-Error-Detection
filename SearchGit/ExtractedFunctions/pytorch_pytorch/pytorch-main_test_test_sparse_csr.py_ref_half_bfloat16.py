        def ref_half_bfloat16(c, a, b, alpha=None, beta=None, out=None):
            res = alpha * (a.to_dense().to(torch.float32) @ b.to_dense().to(torch.float32)).to(a.dtype) + beta * c
            if out is not None:
                out.copy_(res)
                return out
            else:
                return res