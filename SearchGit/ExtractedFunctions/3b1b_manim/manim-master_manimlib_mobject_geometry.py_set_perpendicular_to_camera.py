    def set_perpendicular_to_camera(self, camera_frame):
        to_cam = camera_frame.get_implied_camera_location() - self.get_center()
        normal = self.get_unit_normal()
        axis = normalize(self.get_vector())
        # Project to be perpendicular to axis
        trg_normal = to_cam - np.dot(to_cam, axis) * axis
        mat = rotation_between_vectors(normal, trg_normal)
        self.apply_matrix(mat, about_point=self.get_start())
        return self