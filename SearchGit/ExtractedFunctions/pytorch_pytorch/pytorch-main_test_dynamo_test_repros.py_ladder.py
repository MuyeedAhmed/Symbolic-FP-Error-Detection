        def ladder(x):
            trail = x.size(-1)
            assert trail > 2  # noqa: S101
            weights = []
            for s in [trail, trail - 1, trail - 2]:
                weights.append(torch.ones(s, s - 1))

            for w in weights:
                x = x @ w

            weights.reverse()

            for w in weights:
                x = x @ w.t()

            return x