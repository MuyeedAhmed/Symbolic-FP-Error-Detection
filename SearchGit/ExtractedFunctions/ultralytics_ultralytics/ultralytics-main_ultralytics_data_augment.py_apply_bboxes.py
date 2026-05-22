    def apply_bboxes(self, bboxes: np.ndarray, M: np.ndarray) -> np.ndarray:
        """Apply affine transformation to bounding boxes.

        This function applies an affine transformation to a set of bounding boxes using the provided transformation
        matrix.

        Args:
            bboxes (np.ndarray): Bounding boxes in xyxy format with shape (N, 4), where N is the number of bounding
                boxes.
            M (np.ndarray): Affine transformation matrix with shape (3, 3).

        Returns:
            (np.ndarray): Transformed bounding boxes in xyxy format with shape (N, 4).

        Examples:
            >>> rp = RandomPerspective()
            >>> bboxes = np.array([[10, 10, 20, 20], [30, 30, 40, 40]], dtype=np.float32)
            >>> M = np.eye(3, dtype=np.float32)
            >>> transformed_bboxes = rp.apply_bboxes(bboxes, M)
        """
        n = len(bboxes)
        if n == 0:
            return bboxes

        xy = np.ones((n * 4, 3), dtype=bboxes.dtype)
        xy[:, :2] = bboxes[:, [0, 1, 2, 3, 0, 3, 2, 1]].reshape(n * 4, 2)  # x1y1, x2y2, x1y2, x2y1
        xy = xy @ M.T  # transform
        xy = (xy[:, :2] / xy[:, 2:3] if self.perspective else xy[:, :2]).reshape(n, 8)  # perspective rescale or affine

        # Create new boxes
        x = xy[:, [0, 2, 4, 6]]
        y = xy[:, [1, 3, 5, 7]]
        return np.concatenate((x.min(1), y.min(1), x.max(1), y.max(1)), dtype=bboxes.dtype).reshape(4, n).T