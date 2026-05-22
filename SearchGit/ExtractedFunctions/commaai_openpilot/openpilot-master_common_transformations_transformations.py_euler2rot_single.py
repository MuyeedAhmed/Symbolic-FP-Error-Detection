def euler2rot_single(euler):
  """
  Convert Euler angles (roll, pitch, yaw) to a 3x3 rotation matrix.
  Rotation order: Z-Y-X (yaw, pitch, roll).
  """
  phi, theta, psi = euler

  cx, sx = np.cos(phi), np.sin(phi)
  cy, sy = np.cos(theta), np.sin(theta)
  cz, sz = np.cos(psi), np.sin(psi)

  Rx = np.array([[1, 0, 0], [0, cx, -sx], [0, sx, cx]])
  Ry = np.array([[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]])
  Rz = np.array([[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]])

  return Rz @ Ry @ Rx