    def rodrigues(cls, vectors: npt.NDArray[np.float32]) -> npt.NDArray[np.float32]:
        """Perform batch conversion of rotation vectors to rotation matrices

        Parameters
        ----------
        vectors
            The (N, 3, 1) rotation vectors to convert

        Returns
        -------
        The (N, 3, 3) rotation matrices
        """
        vectors = vectors.reshape(-1, 3)
        theta = np.linalg.norm(vectors, axis=1, keepdims=True)
        units = vectors / (theta + 1e-12)

        k = np.zeros((vectors.shape[0], 3, 3), dtype="float32")
        k[:, 0, 1] = -units[:, 2]
        k[:, 0, 2] = units[:, 1]
        k[:, 1, 0] = units[:, 2]
        k[:, 1, 2] = -units[:, 0]
        k[:, 2, 0] = -units[:, 1]
        k[:, 2, 1] = units[:, 0]

        ident = np.eye(3, dtype="float32")
        retval = ident + np.sin(theta)[:, None] * k + (1 - np.cos(theta))[:, None] * (k @ k)
        return retval