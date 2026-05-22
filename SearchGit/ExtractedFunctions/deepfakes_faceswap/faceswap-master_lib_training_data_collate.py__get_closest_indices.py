    def _get_closest_indices(self) -> tuple[npt.NDArray[np.int64], npt.NDArray[np.int64]]:
        """Obtain the closest x number of landmarks from the opposite side

        Returns
        -------
        indices_a
            Array of size (len(landmarks_a), x) closest B landmarks for each A landmarks
        indices_b
            Array of size (len(landmarks_b), x) closest A landmarks for each B landmarks
        """
        a_count = self._landmarks[0].shape[0]
        b_count = self._landmarks[1].shape[0]
        lms_a = self._landmarks[0].reshape(a_count, -1)
        lms_b = self._landmarks[1].reshape(b_count, -1)

        a_sq = (lms_a ** 2).sum(axis=1, keepdims=True)
        b_sq = (lms_b ** 2).sum(axis=1, keepdims=True)
        dist2 = a_sq + b_sq.T - 2.0 * (lms_a @ lms_b.T)
        np.clip(dist2, 0, None, out=dist2)
        matches_a = np.argpartition(dist2, self._num_choices, axis=1)[:, :self._num_choices]
        matches_b = np.argpartition(dist2.T, self._num_choices, axis=1)[:, :self._num_choices]

        logger.debug("[TrainLoader] Closest matches. A: %s, B: %s",
                     format_array(matches_a), format_array(matches_b))
        return matches_a, matches_b