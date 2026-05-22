    def __wrapped__(self, doc: str, query: str, **kwargs) -> float:
        doc, query = _extract_value(doc), _extract_value(query)

        embeddings = self.model.encode(
            [query, doc], normalize_embeddings=True, **kwargs
        )

        return embeddings[0] @ embeddings[1].T