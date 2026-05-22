    def un_rotate(self,
                  indices_angle: npt.NDArray[np.int32],
                  roi: npt.NDArray[np.float32]) -> None:
        """Un-rotate the given bounding boxes for the given angle indices and update in place

        Parameters
        ----------
        indices_angle
            The angle indices that correlate to the angle each roi was rotated to to obtain the
            result
        roi
            Ragged array of (B, N, 4) detected bounding discovered at the corresponding angle
            index
        """
        mask_needs_rotate = indices_angle > 0
        if not np.any(mask_needs_rotate):
            return

        indices_needs_rotate = np.flatnonzero(mask_needs_rotate)
        matrices = self._matrices_inverse[indices_angle[mask_needs_rotate]]

        for pred_idx, mat in zip(indices_needs_rotate, matrices):
            bboxes = roi[pred_idx]
            pts = np.empty((bboxes.shape[0], 4, 2), dtype="float32")
            pts[:, 0] = bboxes[:, [0, 1]]  # lt
            pts[:, 1] = bboxes[:, [2, 1]]  # rt
            pts[:, 2] = bboxes[:, [2, 3]]  # rb
            pts[:, 3] = bboxes[:, [0, 3]]  # lb

            pts = pts @ mat[:, :2].T + mat[:, 2]

            # boxes must align on (x, y) planes
            bboxes[:, 0] = pts[..., 0].min(axis=1)
            bboxes[:, 1] = pts[..., 1].min(axis=1)
            bboxes[:, 2] = pts[..., 0].max(axis=1)
            bboxes[:, 3] = pts[..., 1].max(axis=1)