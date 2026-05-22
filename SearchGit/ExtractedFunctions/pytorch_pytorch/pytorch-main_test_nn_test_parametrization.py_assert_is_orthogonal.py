        def assert_is_orthogonal(X):
            n, k = X.size(-2), X.size(-1)
            if n < k:
                X = X.mT
                n, k = k, n
            Id = torch.eye(k, dtype=X.dtype, device=X.device).expand(
                *(X.size()[:-2]), k, k
            )
            eps = 10 * n * torch.finfo(X.dtype).eps
            torch.testing.assert_close(X.mH @ X, Id, atol=eps, rtol=0.0)