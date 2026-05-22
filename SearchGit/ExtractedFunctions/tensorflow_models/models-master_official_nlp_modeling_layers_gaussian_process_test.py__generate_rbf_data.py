def _generate_rbf_data(x_data, orthogonal=True):
  """Generates high-dim data that is the eigen components of a RBF kernel."""
  k_rbf = exact_gaussian_kernel(x_data, x_data)
  x_orth, x_diag, _ = np.linalg.svd(k_rbf)
  if orthogonal:
    return x_orth
  return np.diag(np.sqrt(x_diag)).dot(x_orth.T)