    def _get_inverse_matrices(cls, matrices: np.ndarray) -> np.ndarray:
        """ Obtain the inverse matrices for the given matrices. If ``None`` is supplied return a
        dummy transformation matrix that performs no action

        Parameters
        ----------
        matrices : :class:`numpy.ndarray`
            The original transform matrices that the inverse needs to be calculated for

        Returns
        -------
        :class:`numpy.ndarray`
            The inverse transformation matrices
        """
        if not np.any(matrices):
            return np.array([[[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]], dtype=np.float32)

        identity = np.array([[[0., 0., 1.]]], dtype=np.float32)
        mat = np.concatenate([matrices, np.repeat(identity, matrices.shape[0], axis=0)], axis=1)
        retval = np.linalg.inv(mat)
        logger.trace("matrix: %s, inverse: %s", mat, retval)  # type:ignore[attr-defined]
        return retval