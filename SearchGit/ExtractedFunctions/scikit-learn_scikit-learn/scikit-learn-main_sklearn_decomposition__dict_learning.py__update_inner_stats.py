    def _update_inner_stats(self, X, code, batch_size, step):
        """Update the inner stats inplace."""
        if step < batch_size - 1:
            theta = (step + 1) * batch_size
        else:
            theta = batch_size**2 + step + 1 - batch_size
        beta = (theta + 1 - batch_size) / (theta + 1)

        self._A *= beta
        self._A += code.T @ code / batch_size
        self._B *= beta
        self._B += X.T @ code / batch_size