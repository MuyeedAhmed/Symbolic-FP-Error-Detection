def learn_whitening(descriptors, qidxs, pidxs):
  """Learning the post-processing of fine-tuned descriptor vectors.

  This method of whitening learning leverages the provided labeled data and
  uses linear discriminant projections. The projection is decomposed into two
  parts: whitening and rotation. The whitening part is the inverse of the
  square-root of the intraclass (matching pairs) covariance matrix. The
  rotation part is the PCA of the interclass (non-matching pairs) covariance
  matrix in the whitened space. The described approach acts as a
  post-processing step, equivalently, once the fine-tuning of the CNN is
  finished. For more information about the method refer to the section 3.4
  of https://arxiv.org/pdf/1711.02512.pdf.

  Args:
    descriptors: [N, D] NumPy array of L2-normalized descriptors.
    qidxs: List of query indexes.
    pidxs: List of positive pairs indexes.

  Returns:
    mean_descriptor_vector: [N, 1] NumPy array, mean descriptor vector.
    projection: [N, N] NumPy array, whitening projection matrix.
  """
  # Calculating the mean descriptor vector, which is used to perform centering.
  mean_descriptor_vector = descriptors[:, qidxs].mean(axis=1, keepdims=True)
  # Interclass (matching pairs) difference.
  interclass_difference = descriptors[:, qidxs] - descriptors[:, pidxs]
  covariance_matrix = (
      np.dot(interclass_difference, interclass_difference.T) /
      interclass_difference.shape[1])

  # Whitening part.
  projection = np.linalg.inv(cholesky(covariance_matrix))

  projected_descriptors = np.dot(projection,
                                 descriptors - mean_descriptor_vector)
  non_matching_covariance_matrix = np.dot(projected_descriptors,
                                          projected_descriptors.T)
  eigval, eigvec = np.linalg.eig(non_matching_covariance_matrix)
  order = eigval.argsort()[::-1]
  eigvec = eigvec[:, order]

  # Rotational part.
  projection = np.dot(eigvec.T, projection)
  return mean_descriptor_vector, projection