  def _ned_from_calib(self, orientation: Measurement):
    ned_from_device = rot_from_euler(orientation.xyz)
    ned_from_calib = ned_from_device @ self.calib_from_device.T
    ned_from_calib_euler_meas = Measurement(euler_from_rot(ned_from_calib), np.full(3, np.nan))
    return ned_from_calib_euler_meas