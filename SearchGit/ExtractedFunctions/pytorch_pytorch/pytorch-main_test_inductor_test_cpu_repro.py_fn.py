        def fn(inp: torch.Tensor):
            I0 = torch.eye(n, dtype=inp.dtype, device=inp.device)
            I = I0.unsqueeze(0).expand(inp.shape[0], n, n).contiguous()
            hermitian = I + 0.5 * (inp @ inp.mH)
            chol = torch.linalg.cholesky(hermitian, upper=True)
            return chol.abs().sum()