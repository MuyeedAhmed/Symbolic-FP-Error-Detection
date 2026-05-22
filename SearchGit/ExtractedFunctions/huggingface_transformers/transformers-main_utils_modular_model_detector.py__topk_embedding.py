    def _topk_embedding(
        self,
        query_embedding_row: np.ndarray,
        base_embeddings: np.ndarray,
        identifier_map: dict[int, str],
        self_model_normalized: str,
        self_name: str,
        k: int,
    ) -> list[tuple[str, float]]:
        similarities = query_embedding_row @ base_embeddings.T
        indices = np.argpartition(-similarities, k + 32)[: k + 32]
        indices = indices[np.argsort(-similarities[indices])]
        output = []
        for match_id in indices:
            identifier = identifier_map[int(match_id)]
            parent_relative_path, match_name = identifier.split(":", 1)
            parent_model = Path(parent_relative_path).parts[0]
            if match_name == self_name:
                continue
            if self_model_normalized and _normalize(parent_model) == self_model_normalized:
                continue
            output.append((identifier, float(similarities[match_id])))
            if len(output) >= k:
                break
        return output