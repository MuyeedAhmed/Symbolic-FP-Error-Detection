        def fct_no_grad_bhh_Whh(x: torch.Tensor, y: torch.Tensor):
            c_new = y @ W_ih + b_ih + x

            h_new = c_new + 1
            with torch.no_grad():
                h_new_no_grad = torch.tanh(x @ W_hh + b_hh)
            h_new2 = h_new + h_new_no_grad

            return c_new, h_new2