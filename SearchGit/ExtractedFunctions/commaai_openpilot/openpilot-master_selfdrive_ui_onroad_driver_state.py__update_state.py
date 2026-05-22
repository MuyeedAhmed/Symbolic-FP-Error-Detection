  def _update_state(self):
    """Update the driver monitoring state based on model data"""
    sm = ui_state.sm
    if not self.is_visible:
      return

    # Get monitoring state
    dm_state = sm["driverMonitoringState"]
    self.is_active = dm_state.activePolicy == log.DriverMonitoringState.MonitoringPolicy.vision
    self.is_rhd = dm_state.isRHD

    # Update fade state (smoother transition between active/inactive)
    fade_target = 0.0 if self.is_active else 0.5
    self.dm_fade_state = np.clip(self.dm_fade_state + 0.2 * (fade_target - self.dm_fade_state), 0.0, 1.0)

    # Get driver orientation data from appropriate camera
    driverstate = sm["driverStateV2"]
    driver_data = driverstate.rightDriverData if self.is_rhd else driverstate.leftDriverData
    driver_orient = driver_data.faceOrientation

    # Update pose values with scaling and smoothing
    driver_orient = np.array(driver_orient)
    scales = np.where(driver_orient < 0, SCALES_NEG, SCALES_POS)
    v_this = driver_orient * scales
    self.driver_pose_diff = np.abs(self.driver_pose_vals - v_this)
    self.driver_pose_vals = 0.8 * v_this + 0.2 * self.driver_pose_vals  # Smooth changes

    # Apply fade to rotation and compute sin/cos
    rotation_amount = self.driver_pose_vals * (1.0 - self.dm_fade_state)
    self.driver_pose_sins = np.sin(rotation_amount)
    self.driver_pose_coss = np.cos(rotation_amount)

    # Create rotation matrix for 3D face model
    sin_y, sin_x, sin_z = self.driver_pose_sins
    cos_y, cos_x, cos_z = self.driver_pose_coss
    r_xyz = np.array(
      [
        [cos_x * cos_z, cos_x * sin_z, -sin_x],
        [-sin_y * sin_x * cos_z - cos_y * sin_z, -sin_y * sin_x * sin_z + cos_y * cos_z, -sin_y * cos_x],
        [cos_y * sin_x * cos_z - sin_y * sin_z, cos_y * sin_x * sin_z + sin_y * cos_z, cos_y * cos_x],
      ]
    )

    # Transform face keypoints using vectorized matrix multiplication
    self.face_kpts_draw = DEFAULT_FACE_KPTS_3D @ r_xyz.T
    self.face_kpts_draw[:, 2] = self.face_kpts_draw[:, 2] * (1.0 - self.dm_fade_state) + 8 * self.dm_fade_state

    # Pre-calculate the transformed keypoints
    kp_depth = (self.face_kpts_draw[:, 2] - 8) / 120.0 + 1.0
    self.face_keypoints_transformed = self.face_kpts_draw[:, :2] * kp_depth[:, None]

    # Pre-calculate all drawing elements
    self._pre_calculate_drawing_elements()