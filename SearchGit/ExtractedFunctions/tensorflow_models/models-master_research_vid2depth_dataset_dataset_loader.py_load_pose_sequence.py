  def load_pose_sequence(self, frames, target_index):
    """Returns a sequence of pose vectors for frames around the target frame."""
    target_drive, _, target_frame_id = frames[target_index].split(' ')
    target_pose = self.load_pose_raw(target_drive, target_frame_id)
    start_index, end_index = get_seq_start_end(target_frame_id, self.seq_length)
    pose_seq = []
    for index in range(start_index, end_index + 1):
      if index == target_frame_id:
        continue
      drive, _, frame_id = frames[index].split(' ')
      pose = self.load_pose_raw(drive, frame_id)
      # From target to index.
      pose = np.dot(np.linalg.inv(pose), target_pose)
      pose_seq.append(pose)
    return pose_seq