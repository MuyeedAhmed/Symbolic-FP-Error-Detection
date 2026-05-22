  def intrinsics_inv(self):
    # aka 'K_inv' aka view_frame_from_camera_frame
    return np.linalg.inv(self.intrinsics)