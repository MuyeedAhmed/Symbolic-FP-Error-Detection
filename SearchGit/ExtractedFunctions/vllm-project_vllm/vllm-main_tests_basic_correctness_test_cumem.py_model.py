    def model(x):
        out = x @ weight
        cache[: out.size(0)].copy_(out)
        return out + 1