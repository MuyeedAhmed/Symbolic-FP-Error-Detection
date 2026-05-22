    def _rgb_xyz_rgb(self, image: torch.Tensor, mapping: torch.Tensor) -> torch.Tensor:
        """RGB to XYZ or XYZ to RGB conversion.

        Notes
        -----
        The conversion in both directions is the same, but the mapping matrix for XYZ to RGB is
        the inverse of RGB to XYZ.

        References
        ----------
        https://www.image-engineering.de/library/technotes/958-how-to-convert-between-srgb-and-ciexyz

        Parameters
        ----------
        mapping
            The mapping matrix to perform either the XYZ to RGB or RGB to XYZ color space
            conversion

        image
            The image tensor in RGB format

        Returns
        -------
        The image tensor in XYZ format
        """
        dim = image.shape
        image = image.reshape(dim[0], dim[1], dim[2] * dim[3])
        converted = mapping @ image
        return converted.view(dim)