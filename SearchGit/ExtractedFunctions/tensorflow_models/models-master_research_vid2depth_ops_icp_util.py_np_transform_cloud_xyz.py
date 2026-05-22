def np_transform_cloud_xyz(cloud, transform):
  num_points = cloud.shape[0]
  ones = np.ones(shape=[num_points, 1], dtype=np.float32)
  hom_cloud = np.concatenate((cloud, ones), axis=1)
  hom_cloud_t = np.transpose(hom_cloud)
  mat = np_get_transformation_matrix(transform)
  transformed_cloud = np.dot(mat, hom_cloud_t)
  transformed_cloud = np.transpose(transformed_cloud)
  transformed_cloud = transformed_cloud[:, :3]
  return transformed_cloud