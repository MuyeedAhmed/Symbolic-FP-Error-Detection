    def get_images(self,  # pylint:disable=too-many-locals
                   matrices: npt.NDArray[np.float32],
                   feeds: int) -> npt.NDArray[np.float32]:
        """Obtain the sub-crops from the main image patches based on the roi stored in the batch
        and populate them to the batch's data attribute

        Parameters
        ----------
        matrices
            The matrices that define the crops to extract from the expanded patch in shape
            (N x total_feeds, 3, 3)
        feeds
            The number of feeds that are to be made through the model for this batch

        Returns
        -------
        The aligned images that are to be used for 2nd pass re-align
        """
        mats = matrices.reshape(-1, feeds, 3, 3)
        all_offsets = np.rint(mats[..., :2, 2]).astype("int32")
        all_scales = mats[..., 0, 0]  # Always same x/y scaling, always aligned
        all_interpolations = np.where(all_scales < 1.0, cv2.INTER_CUBIC, cv2.INTER_AREA)
        all_dims = np.rint(self._size / all_scales).astype(np.int32)  # Always square

        size = (self._size, self._size)
        retval = np.empty((*mats.shape[:2], *size, 3), dtype=self._images.dtype)

        for batch_id, (offsets, scales, interpolations, dims) in enumerate(zip(all_offsets,
                                                                               all_scales,
                                                                               all_interpolations,
                                                                               all_dims)):
            img = self._images[batch_id]
            for feed_id, offset in enumerate(offsets):
                scale = scales[feed_id]
                interpolation = interpolations[feed_id]
                src_dim = dims[feed_id]
                crop = img[offset[1]:offset[1] + src_dim, offset[0]:offset[0] + src_dim]
                if scale != 1.:
                    crop = cv2.resize(crop, size, interpolation=interpolation)
                retval[batch_id, feed_id] = crop

        # Add the adjusted matrices to :attr:`_matrices` for warping back to frame downstream
        base_mats = self._matrices.reshape(self._matrices.shape[0], -1, 3, 3)
        base_mats = base_mats @ mats @ np.diag([self._size, self._size, 1]).astype("float32")
        self._matrices = base_mats.reshape(matrices.shape[0], *base_mats.shape[2:])

        return retval.reshape(matrices.shape[0], *retval.shape[2:])