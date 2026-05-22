    def _rank_leaf_pairs(self, leaves: list[_PsiTreeNode]) -> np.ndarray:
        """Rank all leaf pairs by original embedding-space cosine similarity."""
        node_embeddings = np.asarray([leaf.embedding for leaf in leaves], dtype=np.float64)
        node_embeddings = self._normalize_embeddings(node_embeddings)
        similarities = node_embeddings @ node_embeddings.T
        lower = np.tril_indices(len(leaves), -1)
        ordered = np.argsort(similarities[lower], axis=0)[::-1]
        return np.stack([lower[0][ordered], lower[1][ordered]], axis=-1)