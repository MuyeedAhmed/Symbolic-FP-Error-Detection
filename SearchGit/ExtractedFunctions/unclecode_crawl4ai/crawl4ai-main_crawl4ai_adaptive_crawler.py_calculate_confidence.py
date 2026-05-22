    async def calculate_confidence(self, state: CrawlState) -> float:
        """Coverage-based learning score (0–1)."""
        # Guard clauses
        if state.kb_embeddings is None or state.query_embeddings is None:
            return 0.0
        if len(state.kb_embeddings) == 0 or len(state.query_embeddings) == 0:
            return 0.0

        # Prepare L2-normalised arrays
        Q = np.asarray(state.query_embeddings, dtype=np.float32)
        D = np.asarray(state.kb_embeddings, dtype=np.float32)
        Q /= np.linalg.norm(Q, axis=1, keepdims=True) + 1e-8
        D /= np.linalg.norm(D, axis=1, keepdims=True) + 1e-8

        # Best cosine per query
        best = (Q @ D.T).max(axis=1)

        # Mean similarity or hit-rate above tau
        tau = getattr(self.config, 'coverage_tau', None)
        score = float((best >= tau).mean()) if tau is not None else float(best.mean())

        # Store quick metrics
        state.metrics['coverage_score'] = score
        state.metrics['avg_best_similarity'] = float(best.mean())
        state.metrics['median_best_similarity'] = float(np.median(best))

        return score