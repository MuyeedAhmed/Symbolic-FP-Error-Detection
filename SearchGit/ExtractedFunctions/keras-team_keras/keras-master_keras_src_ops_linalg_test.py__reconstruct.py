        def _reconstruct(lu, pivots, m, n):
            lower = np.tril(lu[:, : min(m, n)], -1) + np.eye(m, min(m, n))
            upper = np.triu(lu[: min(m, n)])

            # pivots are defined differently in tensorflow
            # compared to the other backends
            if backend.backend() == "tensorflow":
                p_matrix = np.eye(m)[pivots]
            else:
                p_matrix = _pivot_matrix(pivots, m)
            out = p_matrix @ lower @ upper
            return out