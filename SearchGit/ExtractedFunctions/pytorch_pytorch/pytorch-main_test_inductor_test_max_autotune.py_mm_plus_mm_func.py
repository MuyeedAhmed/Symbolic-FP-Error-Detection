        def mm_plus_mm_func(a1, b1, a2, b2) -> torch.Tensor:
            a1_t = torch.permute(a1, [1, 0]).to(torch.bfloat16)
            b1_dtype = b1.to(torch.bfloat16)
            a2_t = torch.permute(a2, [1, 0]).to(torch.bfloat16)
            b2_dtype = b2.to(torch.bfloat16)
            return (a1_t @ b1_dtype + a2_t @ b2_dtype).to(torch.float32)