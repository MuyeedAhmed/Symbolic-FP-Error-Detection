    def matmul_fuse_pattern_two(inp, w1, w2):
        return (inp @ w1, inp @ w2)