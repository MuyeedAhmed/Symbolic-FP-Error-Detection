    def loss_grad(AB):
        # .astype below is needed to ensure y_true and raw_prediction have the
        # same dtype. With result = np.float64(0) * np.array([1, 2], dtype=np.float32)
        # - in Numpy 2, result.dtype is float64
        # - in Numpy<2, result.dtype is float32
        raw_prediction = -(AB[0] * F + AB[1]).astype(dtype=predictions.dtype)
        l, g = bin_loss.loss_gradient(
            y_true=T,
            raw_prediction=raw_prediction,
            sample_weight=sample_weight,
        )
        loss = l.sum()
        # TODO: Remove casting to np.float64 when minimum supported SciPy is 1.11.2
        # With SciPy >= 1.11.2, the LBFGS implementation will cast to float64
        # https://github.com/scipy/scipy/pull/18825.
        # Here we cast to float64 to support SciPy < 1.11.2
        grad = np.asarray([-g @ F, -g.sum()], dtype=np.float64)
        return loss, grad