    def _get_rgb_xyz_map(cls) -> torch.Tensor:
        """Obtain the mapping and inverse mapping for rgb to xyz color space conversion.

        Returns
        -------
        The mapping and inverse Tensors for rgb to xyz color space conversion
        """
        mapping = np.array([[10135552 / 24577794,  8788810 / 24577794, 4435075 / 24577794],
                            [2613072 / 12288897, 8788810 / 12288897, 887015 / 12288897],
                            [1425312 / 73733382, 8788810 / 73733382, 70074185 / 73733382]])
        inverse = np.linalg.inv(mapping)
        return torch.from_numpy(np.stack([mapping, inverse], axis=0)).float()