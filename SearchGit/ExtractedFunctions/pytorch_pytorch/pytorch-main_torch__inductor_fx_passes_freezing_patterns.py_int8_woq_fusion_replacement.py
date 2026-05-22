    def int8_woq_fusion_replacement(inp, w1, w2, w3, s1, s2, s3):
        cat_w = torch.cat((w1, w2, w3), dim=1)
        cat_s = torch.cat((s1, s2, s3), dim=0)
        mm = (inp @ cat_w).mul(cat_s)
        n1, n2 = w1.size(1), w2.size(1)
        return mm.tensor_split([n1, n1 + n2], dim=-1)