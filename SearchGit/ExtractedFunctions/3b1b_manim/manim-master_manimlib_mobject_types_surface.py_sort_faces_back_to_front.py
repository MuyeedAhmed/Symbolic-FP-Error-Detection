    def sort_faces_back_to_front(self, vect: Vect3 = OUT) -> Self:
        tri_is = self.triangle_indices
        points = self.get_points()

        dots = np.dot(points[tri_is[::3]], vect.T)
        indices = np.argsort(dots)
        for k in range(3):
            tri_is[k::3] = tri_is[k::3][indices]
        return self