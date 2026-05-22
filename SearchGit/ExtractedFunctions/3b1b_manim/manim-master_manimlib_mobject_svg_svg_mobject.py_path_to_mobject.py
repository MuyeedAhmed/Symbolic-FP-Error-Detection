    def path_to_mobject(self, path: se.Path, svg: se.SVG) -> VMobjectFromSVGPath:
        if path.id in svg.objects:
            # If this path reuses a referenced definition (<use>), build the mobject from
            # the original geometry.
            # We apply the transform ourselves so we (1) keep the full precision of the 
            # reference and (2) only store one entry in PATH_TO_POINTS.
            ref_path = svg.objects[path.id]
            mob = VMobjectFromSVGPath(ref_path, **self.path_string_config)
            if 'transform' in path.values:
                matrix = se.Matrix(path.values['transform'])
                rotation = np.array([[matrix.a, matrix.b],
                                     [matrix.c, matrix.d]])
                translation = np.array([[matrix.e, matrix.f]])
                mob.apply_points_function(
                    lambda points: np.concatenate([points[:, :2] @ rotation + translation,
                                                   points[:, [2]]],
                                                  axis=1),
                    about_point=None,
                    about_edge=None,
                    works_on_bounding_box=False)
            return mob
        else:
            return VMobjectFromSVGPath(path, **self.path_string_config)