def ned_euler_from_ecef_single(ecef_init, ecef_pose):
  """
  Convert ECEF Euler angles (roll, pitch, yaw) at a given ECEF origin
  to equivalent NED Euler angles.
  """
  converter = LocalCoord(ecef=ecef_init)

  x0 = np.array([1.0, 0, 0])
  y0 = np.array([0, 1.0, 0])
  z0 = np.array([0, 0, 1.0])

  phi, theta, psi = ecef_pose

  x1 = axis_angle_to_rot(z0, psi) @ x0
  y1 = axis_angle_to_rot(z0, psi) @ y0
  z1 = axis_angle_to_rot(z0, psi) @ z0

  x2 = axis_angle_to_rot(y1, theta) @ x1
  y2 = axis_angle_to_rot(y1, theta) @ y1
  z2 = axis_angle_to_rot(y1, theta) @ z1

  x3 = axis_angle_to_rot(x2, phi) @ x2
  y3 = axis_angle_to_rot(x2, phi) @ y2

  zero = np.array(ecef_init)
  x0 = converter.ned2ecef_single([1, 0, 0]) - zero
  y0 = converter.ned2ecef_single([0, 1, 0]) - zero
  z0 = converter.ned2ecef_single([0, 0, 1]) - zero

  psi_out = np.arctan2(np.dot(x3, y0), np.dot(x3, x0))
  theta_out = np.arctan2(-np.dot(x3, z0), np.sqrt(np.dot(x3, x0)**2 + np.dot(x3, y0)**2))

  y2 = axis_angle_to_rot(z0, psi_out) @ y0
  z2 = axis_angle_to_rot(y2, theta_out) @ z0

  phi_out = np.arctan2(np.dot(y3, z2), np.dot(y3, y2))

  return np.array([phi_out, theta_out, psi_out])