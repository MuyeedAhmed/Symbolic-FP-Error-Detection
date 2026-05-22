    def _compute_similarity(self, original: np.ndarray, generated: np.ndarray) -> float:
        """
        Computes the cosine similarity between two embeddings.
        """
        original = original.reshape(1, -1)
        norm = np.linalg.norm(original) * np.linalg.norm(generated, axis=1)
        return np.dot(generated, original.T).flatten() / norm