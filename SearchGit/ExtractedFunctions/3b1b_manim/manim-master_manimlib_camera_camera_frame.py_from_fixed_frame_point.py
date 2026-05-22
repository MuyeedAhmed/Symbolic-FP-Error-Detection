    def from_fixed_frame_point(self, point: Vect3, relative: bool = False):
        inv_view = self.get_inv_view_matrix()
        point4d = [*point, 0 if relative else 1]
        return np.dot(point4d, inv_view.T)[:3]