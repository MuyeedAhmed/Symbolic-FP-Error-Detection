    def extract_custom_vertices(self, polygon, sharp_angle_thresh=45):
        poly = np.array(polygon)
        n = len(poly)
        res = []
        i = 0
        while i < n:
            previous_point = poly[(i - 1) % n]
            current_point = poly[i]
            next_point = poly[(i + 1) % n]
            vector_1 = previous_point - current_point
            vector_2 = next_point - current_point
            cross_product_value = (vector_1[1] * vector_2[0]) - (vector_1[0] * vector_2[1])
            if cross_product_value < 0:
                angle_cos = np.clip(
                    (vector_1 @ vector_2) / (np.linalg.norm(vector_1) * np.linalg.norm(vector_2)), -1.0, 1.0
                )
                angle = np.degrees(np.arccos(angle_cos))
                if abs(angle - sharp_angle_thresh) < 1:
                    # Calculate the new point based on the direction of two vectors.
                    dir_vec = vector_1 / np.linalg.norm(vector_1) + vector_2 / np.linalg.norm(vector_2)
                    dir_vec = dir_vec / np.linalg.norm(dir_vec)
                    step_size = (np.linalg.norm(vector_1) + np.linalg.norm(vector_2)) / 2
                    new_point = current_point + dir_vec * step_size
                    res.append(tuple(new_point))
                else:
                    res.append(tuple(current_point))
            i += 1
        return res