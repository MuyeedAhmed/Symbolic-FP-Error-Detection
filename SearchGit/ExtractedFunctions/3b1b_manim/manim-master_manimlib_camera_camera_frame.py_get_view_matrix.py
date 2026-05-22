    def get_view_matrix(self, refresh=False):
        """
        Returns a 4x4 for the affine transformation mapping a point
        into the camera's internal coordinate system
        """
        if self._data_has_changed:
            shift = self.id4x4.copy()
            rotation = self.id4x4.copy()

            scale = self.get_scale()
            shift[:3, 3] = -self.get_center()
            rotation[:3, :3] = self.get_inverse_camera_rotation_matrix()
            np.dot(rotation, shift, out=self.view_matrix)
            if scale > 0:
                self.view_matrix[:3, :4] /= scale

        return self.view_matrix