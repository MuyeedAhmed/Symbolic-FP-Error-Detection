  def ComputeSimilarity(self,
                        aggregated_descriptors_1,
                        aggregated_descriptors_2,
                        feature_visual_words_1=None,
                        feature_visual_words_2=None):
    """Computes similarity between aggregated descriptors.

    Args:
      aggregated_descriptors_1: 1-D NumPy array.
      aggregated_descriptors_2: 1-D NumPy array.
      feature_visual_words_1: Used only for ASMK/ASMK* aggregation type. 1-D
        sorted NumPy integer array denoting visual words corresponding to
        `aggregated_descriptors_1`.
      feature_visual_words_2: Used only for ASMK/ASMK* aggregation type. 1-D
        sorted NumPy integer array denoting visual words corresponding to
        `aggregated_descriptors_2`.

    Returns:
      similarity: Float. The larger, the more similar.

    Raises:
      ValueError: If aggregation type is invalid.
    """
    if self._aggregation_type == _VLAD:
      similarity = np.dot(aggregated_descriptors_1, aggregated_descriptors_2)
    elif self._aggregation_type == _ASMK:
      similarity = self._AsmkSimilarity(
          aggregated_descriptors_1,
          aggregated_descriptors_2,
          feature_visual_words_1,
          feature_visual_words_2,
          binarized=False)
    elif self._aggregation_type == _ASMK_STAR:
      similarity = self._AsmkSimilarity(
          aggregated_descriptors_1,
          aggregated_descriptors_2,
          feature_visual_words_1,
          feature_visual_words_2,
          binarized=True)
    else:
      raise ValueError('Invalid aggregation type: %d' % self._aggregation_type)

    return similarity