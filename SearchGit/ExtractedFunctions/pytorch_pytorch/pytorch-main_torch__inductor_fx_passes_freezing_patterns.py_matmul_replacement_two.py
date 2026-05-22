    def matmul_replacement_two(inp, w1, w2):
        weights = (w1, w2)
        cat_t = torch.cat(weights, dim=1)
        mm = inp @ cat_t
        return mm.split([w.size(1) for w in weights], dim=1)