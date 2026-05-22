    def search(
        self,
        query: str,
        vectors: List[float],
        top_k: int = 5,
        filters: Optional[Dict] = None,
    ) -> List[OutputData]:
        """
        Search for similar vectors using cosine similarity.

        Args:
            query (str): Query string (not used in vector search)
            vectors (List[float]): Query vector
            top_k (int): Number of results to return
            filters (Dict, optional): Filters to apply to the search

        Returns:
            List[OutputData]: Search results
        """
        filter_conditions = []
        filter_params = []

        if filters:
            for k, v in filters.items():
                filter_conditions.append("JSON_EXTRACT(payload, %s) = %s")
                filter_params.extend([f"$.{k}", json.dumps(v)])

        filter_clause = "WHERE " + " AND ".join(filter_conditions) if filter_conditions else ""

        # For simplicity, we'll compute cosine similarity in Python
        # In production, you'd want to use MySQL stored procedures or UDFs
        with self._get_cursor() as cur:
            query_sql = f"""
                SELECT id, vector, payload
                FROM `{self.collection_name}`
                {filter_clause}
            """
            cur.execute(query_sql, filter_params)
            results = cur.fetchall()

        # Calculate cosine similarity in Python
        import numpy as np
        query_vec = np.array(vectors)
        scored_results = []

        for row in results:
            vec = np.array(json.loads(row['vector']))
            # Cosine similarity
            similarity = np.dot(query_vec, vec) / (np.linalg.norm(query_vec) * np.linalg.norm(vec))
            distance = 1 - similarity
            scored_results.append((row['id'], distance, row['payload']))

        # Sort by distance and apply limit
        scored_results.sort(key=lambda x: x[1])
        scored_results = scored_results[:top_k]

        return [
            OutputData(id=r[0], score=float(r[1]), payload=json.loads(r[2]) if isinstance(r[2], str) else r[2])
            for r in scored_results
        ]