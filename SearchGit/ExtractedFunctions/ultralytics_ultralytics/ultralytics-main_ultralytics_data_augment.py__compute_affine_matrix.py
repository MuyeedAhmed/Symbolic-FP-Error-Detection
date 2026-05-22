    def _compute_affine_matrix(self, img: np.ndarray, size: tuple[int, int]) -> tuple[np.ndarray, float]:
        """Compute the affine transformation matrix without applying it.

        Args:
            img (np.ndarray): Input image used to determine center and dimensions.
            size (tuple[int, int]): Size of the output image (width, height) used for clipping translation transform.

        Returns:
            (M, scale): 3x3 transformation matrix and scale factor.
        """
        # Center
        C = np.eye(3, dtype=np.float32)
        C[0, 2] = -img.shape[1] / 2  # x translation (pixels)
        C[1, 2] = -img.shape[0] / 2  # y translation (pixels)

        # Perspective
        P = np.eye(3, dtype=np.float32)
        P[2, 0] = random.uniform(-self.perspective, self.perspective)  # x perspective (about y)
        P[2, 1] = random.uniform(-self.perspective, self.perspective)  # y perspective (about x)

        # Rotation and Scale
        R = np.eye(3, dtype=np.float32)
        a = random.uniform(-self.degrees, self.degrees)
        if isinstance(self.scale, (tuple, list)):
            s = random.uniform(self.scale[0], self.scale[1])
        else:
            s = random.uniform(1 - self.scale, 1 + self.scale)
        R[:2] = cv2.getRotationMatrix2D(angle=a, center=(0, 0), scale=s)

        # Shear
        S = np.eye(3, dtype=np.float32)
        S[0, 1] = math.tan(random.uniform(-self.shear, self.shear) * math.pi / 180)  # x shear (deg)
        S[1, 0] = math.tan(random.uniform(-self.shear, self.shear) * math.pi / 180)  # y shear (deg)

        # Translation
        T = np.eye(3, dtype=np.float32)

        T[0, 2] = random.uniform(0.5 - self.translate, 0.5 + self.translate) * size[0]  # x translation (pixels)
        T[1, 2] = random.uniform(0.5 - self.translate, 0.5 + self.translate) * size[1]  # y translation (pixels)

        # Combined rotation matrix
        M = T @ S @ R @ P @ C  # order of operations (right to left) is IMPORTANT
        return M, s