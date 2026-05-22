    def add_arc_to(self, point: Vect3, angle: float, n_components: int | None = None, threshold: float = 1e-3) -> Self:
        self.throw_error_if_no_points()
        if abs(angle) < threshold:
            self.add_line_to(point)
            return self

        # Assign default value for n_components
        if n_components is None:
            n_components = int(np.ceil(8 * abs(angle) / TAU))

        arc_points = quadratic_bezier_points_for_arc(angle, n_components)
        target_vect = point - self.get_end()
        curr_vect = arc_points[-1] - arc_points[0]

        arc_points = arc_points @ rotation_between_vectors(curr_vect, target_vect).T
        arc_points *= get_norm(target_vect) / get_norm(curr_vect)
        arc_points += (self.get_end() - arc_points[0])
        self.append_points(arc_points[1:])
        return self