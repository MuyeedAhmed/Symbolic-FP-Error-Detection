    def __call__(self,
                 matrices: npt.NDArray[np.float32],
                 with_roi: bool = False,
                 size: int = 0
                 ) -> npt.NDArray[np.float32] | tuple[npt.NDArray[np.float32],
                                                      npt.NDArray[np.int32]]:
        """Obtain an array of adjusted norm to frame matrices based on the number of re-feed
        iterations that have been selected and the size of the original ROI.

        Parameters
        ----------
        matrices
            A batch of norm to frame transformation matrices to be randomly adjust for re-feeding
            the model in shape (N, 3, 3)
        with_roi
            ``True`` to also return the adjusted ROIs. Default: ``False``
        size
            The size of the image patch that the matrix creates if it cannot be derived from the
            matrices. Default: `0` (derive from matrices)

        Returns
        -------
        matrices
            The adjusted matrices for taking points from normalized to frame space in shape
            ((Num re_feeds * N) + 1, 3, 3), in frame contiguous order (Na, Nb, Nc, Na1, Nb1,
            Nc1...)
        roi
            The ((Num re_feeds * N) + 1, 4) roi for each adjusted feed. Returned if `with_roi` is
            ``True``
        """
        if self._re_feeds == 0:
            raise NotImplementedError
        size_mat = (np.array([size],
                             dtype="float32") if size != 0 else matrices[:, 0, 0])[:, None, None]

        batch_size = matrices.shape[0]
        d_scales = np.random.uniform(1.0 - self.beta,
                                     1.0 + self.beta,
                                     size=(batch_size, self._re_feeds))
        d_shift = size_mat - np.random.uniform(1.0 - self.beta,
                                               1.0 + self.beta,
                                               size=(batch_size, self._re_feeds, 2)) * size_mat

        mats = np.broadcast_to(matrices[:, None], (batch_size, self.total_feeds, 3, 3)).copy()
        mats[:, 1:, (0, 1), (0, 1)] *= d_scales[:, :, None]
        mats[:, 1:, :2, 2] += d_shift
        mats = mats.reshape(-1, 3, 3)
        if not with_roi:
            logger.trace("re-feed. matrices: %s",  # type: ignore[attr-defined]
                         format_array(mats))
            return mats

        tl_br = np.rint((mats @ self._corners).swapaxes(1, 2))
        roi = np.stack([tl_br[:, 0, 0], tl_br[:, 0, 1], tl_br[:, 1, 0], tl_br[:, 1, 1]],
                       axis=1).astype(np.int32)
        logger.trace("re-feed. matrices: %s, roi: %s",  # type: ignore[attr-defined]
                     format_array(mats), format_array(roi))
        return mats, roi