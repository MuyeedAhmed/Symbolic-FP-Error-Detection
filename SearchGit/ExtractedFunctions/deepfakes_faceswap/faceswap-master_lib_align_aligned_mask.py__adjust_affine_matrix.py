    def _adjust_affine_matrix(self, mask_size: int, affine_matrix: np.ndarray) -> np.ndarray:
        """Adjust the affine matrix for the mask's storage size

        Parameters
        ----------
        mask_size
            The original size of the mask.
        affine_matrix
            The affine matrix to transform the mask at original size to the parent frame.

        Returns
        -------
        affine_matrix
            The affine matrix adjusted for the mask at its stored dimensions.
        """
        zoom = self.stored_size / mask_size
        zoom_mat = np.array([[zoom, 0, 0.], [0, zoom, 0.]])
        adjust_mat = np.dot(zoom_mat, self._matrix_2to3(affine_matrix))
        logger.trace("storage_size: %s, mask_size: %s, zoom: %s, "  # type:ignore[attr-defined]
                     "original matrix: %s, adjusted_matrix: %s", self.stored_size, mask_size, zoom,
                     affine_matrix.shape, adjust_mat.shape)
        return adjust_mat