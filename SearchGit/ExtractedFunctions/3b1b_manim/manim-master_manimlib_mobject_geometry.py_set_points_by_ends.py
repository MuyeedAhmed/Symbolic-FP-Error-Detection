    def set_points_by_ends(
        self,
        start: Vect3,
        end: Vect3,
        buff: float = 0,
        path_arc: float = 0
    ) -> Self:
        vect = end - start
        length = max(get_norm(vect), 1e-8)  # More systematic min?
        unit_vect = normalize(vect)

        # Find the right tip length and thickness
        width, tip_width, tip_length = self.get_key_dimensions(length - buff)

        # Adjust start and end based on buff
        if path_arc == 0:
            start = start + buff * unit_vect
            end = end - buff * unit_vect
        else:
            R = length / 2 / math.sin(path_arc / 2)
            midpoint = 0.5 * (start + end)
            center = midpoint + rotate_vector(0.5 * vect, PI / 2) / math.tan(path_arc / 2)
            sign = 1
            start = center + rotate_vector(start - center, buff / R)
            end = center + rotate_vector(end - center, -buff / R)
            path_arc -= (2 * buff + tip_length) / R
        vect = end - start
        length = get_norm(vect)

        # Find points for the stem, imagining an arrow pointed to the left
        if path_arc == 0:
            points1 = (length - tip_length) * np.array([RIGHT, 0.5 * RIGHT, ORIGIN])
            points1 += width * UP / 2
            points2 = points1[::-1] + width * DOWN
        else:
            # Find arc points
            points1 = quadratic_bezier_points_for_arc(path_arc)
            points2 = np.array(points1[::-1])
            points1 *= (R + width / 2)
            points2 *= (R - width / 2)
            rot_T = rotation_matrix_transpose(PI / 2 - path_arc, OUT)
            for points in points1, points2:
                points[:] = np.dot(points, rot_T)
                points += R * DOWN

        self.set_points(points1)
        # Tip
        self.add_line_to(tip_width * UP / 2)
        self.add_line_to(tip_length * LEFT)
        self.tip_index = len(self.get_points()) - 1
        self.add_line_to(tip_width * DOWN / 2)
        self.add_line_to(points2[0])
        # Close it out
        self.add_subpath(points2)
        self.add_line_to(points1[0])

        # Reposition to match proper start and end
        self.rotate(angle_of_vector(vect) - self.get_angle())
        self.rotate(
            PI / 2 - np.arccos(normalize(vect)[2]),
            axis=rotate_vector(self.get_unit_vector(), -PI / 2),
        )
        self.shift(start - self.get_start())
        return self