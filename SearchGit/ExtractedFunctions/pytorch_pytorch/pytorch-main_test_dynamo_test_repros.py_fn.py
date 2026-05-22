        def fn(
            src: torch.Tensor, count: torch.Tensor
        ) -> tuple[tuple[int, ...], tuple[int, ...]]:
            Q, R = torch.linalg.qr(src)
            rhs = torch.ones(Q.shape[0], 1, device=src.device)
            a = torch.linalg.solve_triangular(R, Q.T @ rhs, upper=True)
            cloned = a.clone(memory_format=torch.preserve_format)
            return a.stride(), cloned.stride()