    def _set_full_frame_meta(self,  # pylint:disable=too-many-locals
                             mask: align.Mask,
                             mask_scale: float) -> None:
        """Sets the meta information for displaying the mask in full frame mode.

        Parameters
        ----------
        mask
            The mask object
        mask_scale
            The scaling factor from the stored mask size to the internal mask size

        Sets the following parameters to :attr:`_meta`:
            - roi_mask: the rectangular ROI box from the full frame that contains the original ROI
            for the full frame mask
            - top_left: The location that the roi_mask should be placed in the display frame
            - affine_matrix: The matrix for transposing the mask to a full frame
            - interpolator: The cv2 interpolation method to use for transposing mask to a
            full frame
            - slices: The (`x`, `y`) slice objects required to extract the mask ROI
            from the full frame
        """
        frame_dims = self._globals.current_frame.display_dims
        scaled_mask_roi = np.rint(mask.original_roi *
                                  self._globals.current_frame.scale).astype("int32")

        # Scale and clip the ROI to fit within display frame boundaries
        clipped_roi = scaled_mask_roi.clip(min=(0, 0), max=frame_dims)

        # Obtain min and max points to get ROI as a rectangle
        min_max = {"min": clipped_roi.min(axis=0), "max": clipped_roi.max(axis=0)}

        # Create a bounding box rectangle ROI
        roi_dims = np.rint((min_max["max"][1] - min_max["min"][1],
                            min_max["max"][0] - min_max["min"][0])).astype("uint16")
        roi_mask = np.zeros(roi_dims, dtype="uint8")[..., None]
        corners = T.cast(T.Sequence[np.ndarray],
                         np.expand_dims(scaled_mask_roi - min_max["min"], axis=0))
        # Block out areas outside of the actual mask ROI polygon
        cv2.fillPoly(roi_mask, corners, 255)
        logger.trace(  # type:ignore[attr-defined]
            "Setting Full Frame mask ROI. shape: %s", roi_mask.shape)

        # obtain the slices for cropping mask from full frame
        xy_slices = (slice(int(round(min_max["min"][1])), int(round(min_max["max"][1]))),
                     slice(int(round(min_max["min"][0])), int(round(min_max["max"][0]))))

        # Adjust affine matrix for internal mask size and display dimensions
        adjustments = (np.array([[mask_scale, 0., 0.], [0., mask_scale, 0.]]),
                       np.array([[1 / self._globals.current_frame.scale, 0., 0.],
                                 [0., 1 / self._globals.current_frame.scale, 0.],
                                 [0., 0., 1.]]))
        in_matrix = np.dot(adjustments[0], mask.affine_matrix)
        affine_matrix = np.dot(in_matrix, adjustments[1])

        # Get the size of the mask roi box in the frame
        side_sizes = (scaled_mask_roi[1][0] - scaled_mask_roi[0][0],
                      scaled_mask_roi[1][1] - scaled_mask_roi[0][1])
        mask_roi_size = (side_sizes[0] ** 2 + side_sizes[1] ** 2) ** 0.5

        self._meta.setdefault("roi_mask", []).append(roi_mask)
        self._meta.setdefault("affine_matrix", []).append(affine_matrix)
        self._meta.setdefault("interpolator", []).append(mask.interpolator)
        self._meta.setdefault("slices", []).append(xy_slices)
        self._meta.setdefault("top_left", []).append(min_max["min"] + self._canvas.offset)
        self._meta.setdefault("mask_roi_size", []).append(mask_roi_size)