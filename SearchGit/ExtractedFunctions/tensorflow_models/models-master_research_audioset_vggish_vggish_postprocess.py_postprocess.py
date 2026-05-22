  def postprocess(self, embeddings_batch):
    """Applies postprocessing to a batch of embeddings.

    Args:
      embeddings_batch: An nparray of shape [batch_size, embedding_size]
        containing output from the embedding layer of VGGish.

    Returns:
      An nparray of the same shape as the input but of type uint8,
      containing the PCA-transformed and quantized version of the input.
    """
    assert len(embeddings_batch.shape) == 2, (
        'Expected 2-d batch, got %r' % (embeddings_batch.shape,))
    assert embeddings_batch.shape[1] == vggish_params.EMBEDDING_SIZE, (
        'Bad batch shape: %r' % (embeddings_batch.shape,))

    # Apply PCA.
    # - Embeddings come in as [batch_size, embedding_size].
    # - Transpose to [embedding_size, batch_size].
    # - Subtract pca_means column vector from each column.
    # - Premultiply by PCA matrix of shape [output_dims, input_dims]
    #   where both are are equal to embedding_size in our case.
    # - Transpose result back to [batch_size, embedding_size].
    pca_applied = np.dot(self._pca_matrix,
                         (embeddings_batch.T - self._pca_means)).T

    # Quantize by:
    # - clipping to [min, max] range
    clipped_embeddings = np.clip(
        pca_applied, vggish_params.QUANTIZE_MIN_VAL,
        vggish_params.QUANTIZE_MAX_VAL)
    # - convert to 8-bit in range [0.0, 255.0]
    quantized_embeddings = (
        (clipped_embeddings - vggish_params.QUANTIZE_MIN_VAL) *
        (255.0 /
         (vggish_params.QUANTIZE_MAX_VAL - vggish_params.QUANTIZE_MIN_VAL)))
    # - cast 8-bit float to uint8
    quantized_embeddings = quantized_embeddings.astype(np.uint8)

    return quantized_embeddings