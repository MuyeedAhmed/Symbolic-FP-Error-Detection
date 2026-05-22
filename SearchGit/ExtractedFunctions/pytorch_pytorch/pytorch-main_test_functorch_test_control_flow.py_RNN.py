            def RNN(x: torch.Tensor, y: torch.Tensor):
                c_new_0 = x[0] + 1
                c_new_1 = x[1] + 1
                h_new = (
                    torch.tanh(c_new_1 + x[0] @ W_hh + b_hh)
                    + y[0] @ W_ih
                    + y[1] @ W_ih
                    + b_ih
                    + x[1]
                )
                return (c_new_0, c_new_1), h_new