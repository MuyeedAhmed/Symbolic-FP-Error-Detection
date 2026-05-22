    def _matvec(self, v):
        v = v.ravel()
        n_features = self.shape[0]
        res = np.empty(n_features, dtype=self.X.dtype)
        res[:-1] = safe_sparse_dot(self.X.T, v, dense_output=True) - (
            self.X_mean * self.sqrt_sw.dot(v)
        )
        res[-1] = np.dot(v, self.sqrt_sw)
        return res