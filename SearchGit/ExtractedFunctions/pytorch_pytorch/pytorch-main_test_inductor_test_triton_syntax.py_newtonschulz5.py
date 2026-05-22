        def newtonschulz5(G, steps: int, eps=1e-7):
            assert len(G.shape) == 2  # noqa: S101
            a, b, c = (3.4445, -4.7750, 2.0315)
            X = G.to(
                torch.bfloat16
                if torch.cuda.is_bf16_supported(including_emulation=False)
                else torch.float16
            )
            X /= X.norm() + eps  # ensure top singular value <= 1
            if G.size(0) > G.size(1):
                X = X.T
            for _ in range(steps):
                A = X @ X.T
                B = b * A + c * A @ A
                X = a * X + B @ X
            if G.size(0) > G.size(1):
                X = X.T
            return X