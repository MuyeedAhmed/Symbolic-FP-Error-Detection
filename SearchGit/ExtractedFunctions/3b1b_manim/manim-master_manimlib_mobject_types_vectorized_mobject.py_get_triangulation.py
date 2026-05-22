    def get_triangulation(self) -> np.ndarray:
        # Figure out how to triangulate the interior to know
        # how to send the points as to the vertex shader.
        # First triangles come directly from the points
        points = self.get_points()

        if len(points) <= 1:
            return np.zeros(0, dtype='i4')

        normal_vector = self.get_unit_normal()

        # Rotate points such that unit normal vector is OUT
        if not np.isclose(normal_vector, OUT).all():
            points = np.dot(points, z_to_vector(normal_vector))

        v01s = points[1::2] - points[0:-1:2]
        v12s = points[2::2] - points[1::2]
        curve_orientations = np.sign(cross2d(v01s, v12s))

        concave_parts = curve_orientations < 0

        # These are the vertices to which we'll apply a polygon triangulation
        indices = np.arange(len(points), dtype=int)
        inner_vert_indices = np.hstack([
            indices[0::2],
            indices[1::2][concave_parts],
        ])
        inner_vert_indices.sort()
        # Even indices correspond to anchors, and `end_indices // 2`
        # shows which anchors are considered end points
        end_indices = self.get_subpath_end_indices()
        counts = np.arange(1, len(inner_vert_indices) + 1)
        rings = counts[inner_vert_indices % 2 == 0][end_indices // 2]

        # Triangulate
        inner_verts = points[inner_vert_indices]
        inner_tri_indices = inner_vert_indices[
            earclip_triangulation(inner_verts, rings)
        ]
        # Remove null triangles, coming from adjascent points
        iti = inner_tri_indices
        null1 = (iti[0::3] + 1 == iti[1::3]) & (iti[0::3] + 2 == iti[2::3])
        null2 = (iti[0::3] - 1 == iti[1::3]) & (iti[0::3] - 2 == iti[2::3])
        inner_tri_indices = iti[~(null1 | null2).repeat(3)]

        ovi = self.get_outer_vert_indices()
        tri_indices = np.hstack([ovi, inner_tri_indices])
        return tri_indices