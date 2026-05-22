def rotate_cov(rot_matrix, cov_in):
  return rot_matrix @ cov_in @ rot_matrix.T