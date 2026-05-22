    def _matmat(self, v):
        n_features = self.shape[0]
        res = np.empty((n_features, v.shape[1]), dtype=self.X.dtype)
        res[:-1] = safe_sparse_dot(self.X.T, v, dense_output=True) - self.X_mean[
            :, None
        ] * self.sqrt_sw.dot(v)
        res[-1] = np.dot(self.sqrt_sw, v)
        return res