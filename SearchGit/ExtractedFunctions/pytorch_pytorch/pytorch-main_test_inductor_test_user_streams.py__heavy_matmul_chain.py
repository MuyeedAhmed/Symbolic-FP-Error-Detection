    def _heavy_matmul_chain(x, w, depth=8):
        """Chain of matmuls to create substantial GPU work (~ms).
        Used to widen the race window between streams so that missing
        synchronization is observable."""
        h = x
        for _ in range(depth):
            h = h @ w
        return h