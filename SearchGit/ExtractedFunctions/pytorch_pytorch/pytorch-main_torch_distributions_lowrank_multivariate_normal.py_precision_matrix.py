    def precision_matrix(self) -> Tensor:
        # We use "Woodbury matrix identity" to take advantage of low rank form::
        #     inv(W @ W.T + D) = inv(D) - inv(D) @ W @ inv(C) @ W.T @ inv(D)
        # where :math:`C` is the capacitance matrix.
        Wt_Dinv = (
            self._unbroadcasted_cov_factor.mT
            / self._unbroadcasted_cov_diag.unsqueeze(-2)
        )
        A = torch.linalg.solve_triangular(self._capacitance_tril, Wt_Dinv, upper=False)
        precision_matrix = (
            torch.diag_embed(self._unbroadcasted_cov_diag.reciprocal()) - A.mT @ A
        )
        return precision_matrix.expand(
            self._batch_shape + self._event_shape + self._event_shape
        )