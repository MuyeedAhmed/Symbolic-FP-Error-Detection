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
        try:
            # Fetch all vectors (in production, you'd want pagination or filtering)
            query_cql = f"""
                SELECT id, vector, payload
                FROM {self.keyspace}.{self.collection_name}
            """
            rows = self.session.execute(query_cql)

            # Calculate cosine similarity in Python
            query_vec = np.array(vectors)
            scored_results = []

            for row in rows:
                if not row.vector:
                    continue

                vec = np.array(row.vector)
                
                # Cosine similarity
                similarity = np.dot(query_vec, vec) / (np.linalg.norm(query_vec) * np.linalg.norm(vec))
                distance = 1 - similarity

                # Apply filters if provided
                if filters:
                    try:
                        payload = json.loads(row.payload) if row.payload else {}
                        match = all(payload.get(k) == v for k, v in filters.items())
                        if not match:
                            continue
                    except json.JSONDecodeError:
                        continue

                scored_results.append((row.id, distance, row.payload))

            # Sort by distance and apply limit
            scored_results.sort(key=lambda x: x[1])
            scored_results = scored_results[:top_k]

            return [
                OutputData(
                    id=r[0],
                    score=float(r[1]),
                    payload=json.loads(r[2]) if r[2] else {}
                )
                for r in scored_results
            ]
        except Exception as e:
            logger.error(f"Search failed: {e}")
            raise