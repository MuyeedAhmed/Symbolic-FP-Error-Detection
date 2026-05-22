    def zeropower_via_newtonschulz5(self, x, steps: int):
        """We apply the Newton-Schulz iteration to compute matrix G.

        We select a quintic iteration that maximizes the slope at zero. This
        approach helps minimize steps, even if the iteration doesn't fully
        converge across the interval. The result isn't exactly UV^T (from the
        SVD of G), but rather an approximation like US'V^T. Despite this
        approximation, model performance remains unaffected compared to using
        the exact UV^T from the SVD.
        """
        shape = ops.shape(x)
        if len(shape) < 2:
            raise ValueError(
                "Expected gradient or momentum to have at least 2 dimensions. "
                f"Received: shape={shape}"
            )

        a, b, c = self.muon_a, self.muon_b, self.muon_c
        if shape[-2] > shape[-1]:
            x = self.transpose_last_axis(x)

        # Ensure spectral norm is at most 1
        x = x / (ops.norm(x, axis=(-2, -1), keepdims=True) + 1e-7)
        # Perform the NS iterations
        for _ in range(steps):
            temp_a = x @ self.transpose_last_axis(x)
            temp_b = b * temp_a + c * temp_a @ temp_a
            x = a * x + temp_b @ x

        if shape[-2] > shape[-1]:
            x = self.transpose_last_axis(x)
        return x