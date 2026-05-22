        def _pivot_matrix(pivots, n):
            p_matrix = np.eye(n)
            for i, p in enumerate(pivots):
                identity = np.eye(n, n)
                q = identity[i, :].copy()
                identity[i, :] = identity[p, :]
                identity[p, :] = q
                p_matrix = np.dot(p_matrix, identity)
            return p_matrix