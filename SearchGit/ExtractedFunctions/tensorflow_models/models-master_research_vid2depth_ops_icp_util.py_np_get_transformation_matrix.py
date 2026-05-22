def np_get_transformation_matrix(transform):
  """Converts [tx, ty, tz, rx, ry, rz] to a transform matrix."""
  rx = transform[3]
  ry = transform[4]
  rz = transform[5]

  rz = np.clip(rz, -np.pi, np.pi)
  ry = np.clip(ry, -np.pi, np.pi)
  rx = np.clip(rx, -np.pi, np.pi)

  cos_rx = np.cos(rx)
  sin_rx = np.sin(rx)
  rotx_1 = np.stack([1.0, 0.0, 0.0])
  rotx_2 = np.stack([0.0, cos_rx, -sin_rx])
  rotx_3 = np.stack([0.0, sin_rx, cos_rx])
  xmat = np.stack([rotx_1, rotx_2, rotx_3])

  cos_ry = np.cos(ry)
  sin_ry = np.sin(ry)
  roty_1 = np.stack([cos_ry, 0.0, sin_ry])
  roty_2 = np.stack([0.0, 1.0, 0.0])
  roty_3 = np.stack([-sin_ry, 0.0, cos_ry])
  ymat = np.stack([roty_1, roty_2, roty_3])

  cos_rz = np.cos(rz)
  sin_rz = np.sin(rz)
  rotz_1 = np.stack([cos_rz, -sin_rz, 0.0])
  rotz_2 = np.stack([sin_rz, cos_rz, 0.0])
  rotz_3 = np.stack([0.0, 0.0, 1.0])
  zmat = np.stack([rotz_1, rotz_2, rotz_3])

  rotate = np.dot(np.dot(xmat, ymat), zmat)

  translate = transform[:3]
  mat = np.concatenate((rotate, np.expand_dims(translate, 1)), axis=1)

  hom_filler = np.array([[0.0, 0.0, 0.0, 1.0]], dtype=np.float32)
  mat = np.concatenate((mat, hom_filler), axis=0)
  return mat