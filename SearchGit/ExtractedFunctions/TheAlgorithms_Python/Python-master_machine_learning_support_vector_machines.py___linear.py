    def __linear(self, vector1: ndarray, vector2: ndarray) -> float:
        """Linear kernel (as if no kernel used at all)"""
        return np.dot(vector1, vector2)