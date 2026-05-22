        def f1(x: torch.Tensor, y: torch.Tensor):
            c_new = y @ W_1 + b_1
            h_new = torch.tanh(c_new + x)
            return c_new, h_new