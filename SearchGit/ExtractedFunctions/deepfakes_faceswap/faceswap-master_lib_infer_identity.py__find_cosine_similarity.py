    def _find_cosine_similarity(cls,
                                source: npt.NDArray[np.float32],
                                batch: npt.NDArray[np.float32]) -> npt.NDArray[np.float64]:
        """Find the cosine similarity between a source face identity and a test face identity

        Parameters
        ---------
        source
            The identity encoding for the source face identities
        batch
            A batch of face identities to test against the sources

        Returns
        -------
        The cosine similarity between the face identities and the source identities
        """
        s_norms = source / np.linalg.norm(source, axis=1, keepdims=True)
        t_norms = batch / np.linalg.norm(batch, axis=1, keepdims=True)
        retval = t_norms @ s_norms.T
        return retval