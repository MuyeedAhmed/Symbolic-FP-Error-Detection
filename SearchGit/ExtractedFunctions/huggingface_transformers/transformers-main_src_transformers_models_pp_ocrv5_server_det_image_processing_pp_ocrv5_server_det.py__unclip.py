    def _unclip(self, contour_box, unclip_ratio):
        """
        Expands (dilates) a detected text bounding box to recover the full text region.

        Args:
            contour_box (np.ndarray): Input contour of shape (N, 2), where N is the number of points.
            unclip_ratio (float): Expansion ratio, typically greater than 1.0.

        Returns:
            np.ndarray: Expanded contour of shape (M, 2).
        """
        # --- 1. Parameter calculation ---
        polygon = contour_box.reshape(-1, 2).astype(np.float32)
        perimeter = cv2.arcLength(polygon, True)
        area = cv2.contourArea(polygon)
        offset_distance = area * unclip_ratio / perimeter

        # --- 2. Determine polygon orientation and edge normals ---
        x, y = polygon[:, 0], polygon[:, 1]
        is_counter_clockwise = (x @ np.roll(y, -1) - y @ np.roll(x, -1)) > 0.0

        edges = np.roll(polygon, -1, axis=0) - polygon
        edge_lengths = np.linalg.norm(edges, axis=1, keepdims=True)
        edge_directions = edges / np.maximum(edge_lengths, 1e-6)

        if is_counter_clockwise:
            normals = np.stack([edge_directions[:, 1], -edge_directions[:, 0]], axis=1)
        else:
            normals = np.stack([-edge_directions[:, 1], edge_directions[:, 0]], axis=1)

        # --- 3. Calculate new vertices from intersecting shifted edge lines ---
        shifted_points = polygon + offset_distance * normals

        prev_shifted_points = np.roll(shifted_points, 1, axis=0)
        prev_edge_directions = np.roll(edge_directions, 1, axis=0)

        cross_product = (
            prev_edge_directions[:, 0] * edge_directions[:, 1] - prev_edge_directions[:, 1] * edge_directions[:, 0]
        )

        is_parallel_mask = np.abs(cross_product) < 1e-6
        cross_product_safe = np.where(is_parallel_mask, 1.0, cross_product)

        vec_to_current = shifted_points - prev_shifted_points
        intersection_param = (
            vec_to_current[:, 0] * edge_directions[:, 1] - vec_to_current[:, 1] * edge_directions[:, 0]
        ) / cross_product_safe

        new_vertices = prev_shifted_points + prev_edge_directions * intersection_param[:, None]

        # --- 4. Handle near-parallel adjacent edges with a fallback ---
        if np.any(is_parallel_mask):
            prev_normals = np.roll(normals, 1, axis=0)
            fallback_points = polygon + 0.5 * offset_distance * (prev_normals + normals)
            new_vertices[is_parallel_mask] = fallback_points[is_parallel_mask]

        return np.array([new_vertices.astype(np.float32)])