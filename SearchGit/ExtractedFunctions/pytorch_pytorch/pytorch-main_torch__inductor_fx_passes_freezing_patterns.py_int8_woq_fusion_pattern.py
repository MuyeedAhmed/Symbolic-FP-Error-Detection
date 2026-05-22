    def int8_woq_fusion_pattern(inp, w1, w2, w3, s1, s2, s3):
        return ((inp @ w1) * s1, (inp @ w2) * s2, (inp @ w3) * s3)