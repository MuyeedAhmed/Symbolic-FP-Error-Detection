    def covariance_matrix(self) -> Tensor:
        return (
            self._unbroadcasted_scale_tril
            @ self._unbroadcasted_scale_tril.transpose(-2, -1)
        ).expand(self._batch_shape + self._event_shape)