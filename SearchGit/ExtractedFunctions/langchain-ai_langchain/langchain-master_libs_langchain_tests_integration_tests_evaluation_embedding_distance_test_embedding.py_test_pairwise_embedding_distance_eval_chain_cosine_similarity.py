def test_pairwise_embedding_distance_eval_chain_cosine_similarity(
    pairwise_embedding_distance_eval_chain: PairwiseEmbeddingDistanceEvalChain,
    vectors: tuple[np.ndarray, np.ndarray],
) -> None:
    """Test the cosine similarity."""
    pairwise_embedding_distance_eval_chain.distance_metric = EmbeddingDistance.COSINE
    result = pairwise_embedding_distance_eval_chain._compute_score(np.array(vectors))
    expected = 1.0 - np.dot(vectors[0], vectors[1]) / (
        np.linalg.norm(vectors[0]) * np.linalg.norm(vectors[1])
    )
    assert np.isclose(result, expected)