  def _transform_calib_from_device(self, meas: Measurement):
    new_xyz = self.calib_from_device @ meas.xyz
    new_xyz_std = rotate_std(self.calib_from_device, meas.xyz_std)
    return Measurement(new_xyz, new_xyz_std)