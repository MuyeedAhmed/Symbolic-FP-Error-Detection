        def mm_func(a, b) -> torch.Tensor:
            a_t = torch.permute(a, [1, 0]).to(torch.bfloat16)
            b_dtype = b.to(torch.bfloat16)
            # Add .to() to make sure that mm could be potentially padded
            # Strides for output are not padded
            return (a_t @ b_dtype).to(torch.float32)