        def second_chain_fct(scan_fct, inp, **kwargs):
            W = torch.ones(2, 5, device=device)
            return inp @ W