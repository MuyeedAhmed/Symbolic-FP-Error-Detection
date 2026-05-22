    def test_affine_2d_rotateRandom(self, device):
        # scipy before 1.0.0 do not support homogeneous coordinate
        # scipy.ndimage.affine_transform, so we need to skip.
        for angle_rad, input_size2d, output_size2d in \
                itertools.product(angle_rad_(), input_size2d_(), output_size2d_()):

            input_size = input_size2d
            input_ary = np.array(np.random.random(input_size), dtype=np.float32).round(3)
            output_size = output_size2d

            input_ary[0, 0, 0, 0] = 2
            input_ary[0, 0, 0, -1] = 4
            input_ary[0, 0, -1, 0] = 6
            input_ary[0, 0, -1, -1] = 8

            transform_tensor, transform_ary, grid_ary = \
                _buildEquivalentAffineTransforms2d(device, input_size, output_size, angle_rad)

            scipy_ary = torch.from_numpy(scipy.ndimage.affine_transform(
                input_ary[0, 0],
                transform_ary,
                output_shape=output_size[2:],
                order=1,
                mode='nearest',
                prefilter=False))

            affine_tensor = torch.nn.functional.affine_grid(
                transform_tensor,
                torch.Size(output_size),
                align_corners=True
            )

            gridsample_ary = torch.nn.functional.grid_sample(
                torch.tensor(input_ary, device=device).to(device),
                affine_tensor,
                padding_mode='border',
                align_corners=True
            ).to('cpu')

            affine_tensor = affine_tensor.to('cpu')

            for r in range(affine_tensor.size(1)):
                for c in range(affine_tensor.size(2)):
                    grid_out = np.dot(grid_ary, [r, c, 1])
                    self.assertEqual(affine_tensor[0, r, c], grid_out[:2], exact_dtype=False)

            self.assertEqual(scipy_ary, gridsample_ary.reshape_as(scipy_ary))