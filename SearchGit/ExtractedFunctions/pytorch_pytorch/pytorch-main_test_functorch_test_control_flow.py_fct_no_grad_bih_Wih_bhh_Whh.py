        def fct_no_grad_bih_Wih_bhh_Whh(x: torch.Tensor, y: torch.Tensor):
            c_new = x + y
            h_new = c_new + 1
            with torch.no_grad():
                c_new_no_grad = y @ W_ih + b_ih
                h_new_no_grad = torch.tanh(x @ W_hh + b_hh)
            c_new2 = c_new + c_new_no_grad
            h_new2 = h_new + h_new_no_grad
            return c_new2, h_new2