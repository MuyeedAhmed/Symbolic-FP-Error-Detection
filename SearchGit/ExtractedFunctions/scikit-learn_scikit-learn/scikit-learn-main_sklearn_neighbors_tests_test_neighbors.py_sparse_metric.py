    def sparse_metric(x, y):  # Metric accepting sparse matrix input (only)
        assert issparse(x) and issparse(y)
        return x.dot(y.T).toarray().item()