    def _get_transformed_coordinates(
        self, x_coords, y_coords, transformation_matrix
    ):
        backend = self.backend

        batch_size = backend.shape(transformation_matrix)[0]

        homogeneous_transform = backend.numpy.concatenate(
            [transformation_matrix, backend.numpy.ones((batch_size, 1, 1))],
            axis=-1,
        )
        homogeneous_transform = backend.numpy.reshape(
            homogeneous_transform, (batch_size, 3, 3)
        )

        inverse_transform = backend.linalg.inv(homogeneous_transform)

        ones_column = backend.numpy.ones_like(x_coords)
        homogeneous_coords = backend.numpy.concatenate(
            [x_coords, y_coords, ones_column], axis=-1
        )

        homogeneous_coords = backend.numpy.moveaxis(homogeneous_coords, -1, -2)
        transformed_coords = backend.numpy.matmul(
            inverse_transform, homogeneous_coords
        )
        transformed_coords = backend.numpy.moveaxis(transformed_coords, -1, -2)

        x_transformed = transformed_coords[..., 0] / transformed_coords[..., 2]
        y_transformed = transformed_coords[..., 1] / transformed_coords[..., 2]

        return x_transformed, y_transformed