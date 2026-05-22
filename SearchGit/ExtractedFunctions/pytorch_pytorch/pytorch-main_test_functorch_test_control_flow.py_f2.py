        def f2(x: torch.Tensor, y: torch.Tensor):
            c_new = y @ W_2 + b_2 * b_1 + y @ W_1
            h_new = torch.tanh(c_new + x)
            return c_new, h_new