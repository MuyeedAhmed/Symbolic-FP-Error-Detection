def apply_whitening(descriptors,
                    mean_descriptor_vector,
                    projection,
                    output_dim=None):
  """Applies the whitening to the descriptors as a post-processing step.

  Args:
    descriptors: [N, D] NumPy array of L2-normalized descriptors to be
      post-processed.
    mean_descriptor_vector: Mean descriptor vector.
    projection: Whitening projection matrix.
    output_dim: Integer, parameter for the dimensionality reduction. If
      `output_dim` is None, the dimensionality reduction is not performed.

  Returns:
    descriptors_whitened: [N, output_dim] NumPy array of L2-normalized
      descriptors `descriptors` after whitening application.
  """
  eps = 1e-6
  if output_dim is None:
    output_dim = projection.shape[0]

  descriptors = np.dot(projection[:output_dim, :],
                       descriptors - mean_descriptor_vector)
  descriptors_whitened = descriptors / (
      np.linalg.norm(descriptors, ord=2, axis=0, keepdims=True) + eps)
  return descriptors_whitened