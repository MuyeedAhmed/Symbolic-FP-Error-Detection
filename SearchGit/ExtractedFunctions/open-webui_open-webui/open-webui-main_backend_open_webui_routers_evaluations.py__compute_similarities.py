def _compute_similarities(feedbacks: list[LeaderboardFeedbackData], query: str) -> dict:
    """
    Compute how relevant each feedback is to a search query.

    Uses embeddings to find semantic similarity between the query and
    each feedback's tags. Higher similarity means the feedback is more
    relevant to what the user searched for.

    This is used to weight Elo calculations - feedbacks matching the
    query have more influence on the final rankings.

    Returns: {feedback_id: similarity_score (0-1)}
    """
    import numpy as np

    embedding_model = _get_embedding_model()
    if not embedding_model:
        return {}

    all_tags = list({tag for feedback in feedbacks if feedback.data for tag in feedback.data.get('tags', [])})
    if not all_tags:
        return {}

    try:
        tag_embeddings = embedding_model.encode(all_tags)
        query_embedding = embedding_model.encode([query])[0]
    except Exception as e:
        log.error(f'Embedding error: {e}')
        return {}

    # Vectorized cosine similarity
    tag_norms = np.linalg.norm(tag_embeddings, axis=1)
    query_norm = np.linalg.norm(query_embedding)
    similarities = np.dot(tag_embeddings, query_embedding) / (tag_norms * query_norm + 1e-9)
    tag_similarity_map = dict(zip(all_tags, similarities.tolist()))

    return {
        feedback.id: max(
            (tag_similarity_map.get(tag, 0) for tag in (feedback.data or {}).get('tags', [])),
            default=0,
        )
        for feedback in feedbacks
    }