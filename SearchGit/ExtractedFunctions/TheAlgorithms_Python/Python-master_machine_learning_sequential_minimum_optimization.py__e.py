    def _e(self, index):
        """
        Two cases:
            1: Sample[index] is non-bound, fetch error from list: _error
            2: sample[index] is bound, use predicted value minus true value: g(xi) - yi
        """
        # get from error data
        if self._is_unbound(index):
            return self._error[index]
        # get by g(xi) - yi
        else:
            gx = np.dot(self.alphas * self.tags, self._K_matrix[:, index]) + self._b
            yi = self.tags[index]
            return gx - yi