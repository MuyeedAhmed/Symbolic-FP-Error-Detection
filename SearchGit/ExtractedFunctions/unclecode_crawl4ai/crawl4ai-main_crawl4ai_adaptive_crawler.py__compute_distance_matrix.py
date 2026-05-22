    def _compute_distance_matrix(self, query_embeddings: Any, kb_embeddings: Any) -> Any:
        """Compute distance matrix using vectorized operations"""
        
        
        if kb_embeddings is None or len(kb_embeddings) == 0:
            return None
            
        # Ensure proper shapes
        if len(query_embeddings.shape) == 1:
            query_embeddings = query_embeddings.reshape(1, -1)
        if len(kb_embeddings.shape) == 1:
            kb_embeddings = kb_embeddings.reshape(1, -1)
            
        # Vectorized cosine distance: 1 - cosine_similarity
        # Normalize vectors
        query_norm = query_embeddings / np.linalg.norm(query_embeddings, axis=1, keepdims=True)
        kb_norm = kb_embeddings / np.linalg.norm(kb_embeddings, axis=1, keepdims=True)
        
        # Compute cosine similarity matrix
        similarity_matrix = np.dot(query_norm, kb_norm.T)
        
        # Convert to distance
        distance_matrix = 1 - similarity_matrix
        
        return distance_matrix