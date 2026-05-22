    def point_to_number(self, point: Vect3 | Vect3Array) -> float | VectN:
        start = self.get_points()[0]
        end = self.get_points()[-1]
        vect = end - start
        proportion = fdiv(
            np.dot(point - start, vect),
            np.dot(end - start, vect),
        )
        return interpolate(self.x_min, self.x_max, proportion)