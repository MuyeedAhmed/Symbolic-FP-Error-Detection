            def is_orthogonal(self, X):
                Id = torch.eye(X.size(0))
                return torch.allclose(X.T @ X, Id, atol=1e-4)