    def adjust_point(self, poly):
        """
        adjust point order.
        """
        point_num = poly.shape[0]
        if point_num == 4:
            len_1 = np.linalg.norm(poly[0] - poly[1])
            len_2 = np.linalg.norm(poly[1] - poly[2])
            len_3 = np.linalg.norm(poly[2] - poly[3])
            len_4 = np.linalg.norm(poly[3] - poly[0])

            if (len_1 + len_3) * 1.5 < (len_2 + len_4):
                poly = poly[[1, 2, 3, 0], :]

        elif point_num > 4:
            vector_1 = poly[0] - poly[1]
            vector_2 = poly[1] - poly[2]
            cos_theta = np.dot(vector_1, vector_2) / (
                np.linalg.norm(vector_1) * np.linalg.norm(vector_2) + 1e-6
            )
            theta = np.arccos(np.round(cos_theta, decimals=4))

            if abs(theta) > (70 / 180 * math.pi):
                index = list(range(1, point_num)) + [0]
                poly = poly[np.array(index), :]
        return poly