    def apply_segments(
        self, segments: np.ndarray, M: np.ndarray, size: tuple[int, int]
    ) -> tuple[np.ndarray, np.ndarray]:
        """Apply affine transformations to segments and generate new bounding boxes.

        This function applies affine transformations to input segments and generates new bounding boxes based on the
        transformed segments. It clips the transformed segments to fit within the new bounding boxes.

        Args:
            segments (np.ndarray): Input segments with shape (N, M, 2), where N is the number of segments and M is the
                number of points in each segment.
            M (np.ndarray): Affine transformation matrix with shape (3, 3).
            size (tuple[int, int]): Size of the output image (width, height) used for clipping the segments.

        Returns:
            bboxes (np.ndarray): New bounding boxes with shape (N, 4) in xyxy format.
            segments (np.ndarray): Transformed and clipped segments with shape (N, M, 2).

        Examples:
            >>> rp = RandomPerspective()
            >>> segments = np.random.rand(10, 500, 2)  # 10 segments with 500 points each
            >>> M = np.eye(3)  # Identity transformation matrix
            >>> new_bboxes, new_segments = rp.apply_segments(segments, M)
        """
        n, num = segments.shape[:2]
        if n == 0:
            return [], segments

        xy = np.ones((n * num, 3), dtype=segments.dtype)
        segments = segments.reshape(-1, 2)
        xy[:, :2] = segments
        xy = xy @ M.T  # transform
        xy = xy[:, :2] / xy[:, 2:3]
        segments = xy.reshape(n, -1, 2)
        bboxes = np.stack([segment2box(xy, size[0], size[1]) for xy in segments], 0)
        segments[..., 0] = segments[..., 0].clip(bboxes[:, 0:1], bboxes[:, 2:3])
        segments[..., 1] = segments[..., 1].clip(bboxes[:, 1:2], bboxes[:, 3:4])
        return bboxes, segments