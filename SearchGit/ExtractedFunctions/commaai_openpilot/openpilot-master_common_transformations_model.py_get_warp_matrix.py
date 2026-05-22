def get_warp_matrix(device_from_calib_euler: np.ndarray, intrinsics: np.ndarray, bigmodel_frame: bool = False) -> np.ndarray:
  calib_from_model = calib_from_sbigmodel if bigmodel_frame else calib_from_medmodel
  device_from_calib = rot_from_euler(device_from_calib_euler)
  camera_from_calib = intrinsics @ view_frame_from_device_frame @ device_from_calib
  warp_matrix: np.ndarray = camera_from_calib @ calib_from_model
  return warp_matrix