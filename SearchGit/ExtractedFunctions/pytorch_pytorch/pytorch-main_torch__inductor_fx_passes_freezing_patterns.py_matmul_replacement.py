    def matmul_replacement(inp, w1, w2, w3):
        weights = (w1, w2, w3)
        cat_t = torch.cat(weights, dim=1)
        mm = inp @ cat_t
        return mm.split([w.size(1) for w in weights], dim=1)