def embedding_distance(tracks: list, detections: list, metric: str = "cosine") -> np.ndarray:
    """Compute distance between tracks and detections based on embeddings.

    Args:
        tracks (list[BOTrack]): List of tracks, where each track contains embedding features.
        detections (list[BOTrack]): List of detections, where each detection contains embedding features.
        metric (str): Metric for distance computation. Supported metrics include 'cosine', 'euclidean', etc.

    Returns:
        (np.ndarray): Cost matrix computed based on embeddings with shape (N, M), where N is the number of tracks and M
            is the number of detections.

    Examples:
        Compute the embedding distance between tracks and detections using cosine metric
        >>> tracks = [BOTrack(...), BOTrack(...)]  # List of track objects with embedding features
        >>> detections = [BOTrack(...), BOTrack(...)]  # List of detection objects with embedding features
        >>> cost_matrix = embedding_distance(tracks, detections, metric="cosine")
    """
    cost_matrix = np.zeros((len(tracks), len(detections)), dtype=np.float32)
    if cost_matrix.size == 0:
        return cost_matrix
    det_features = np.asarray([track.curr_feat for track in detections], dtype=np.float32)
    # for i, track in enumerate(tracks):
    # cost_matrix[i, :] = np.maximum(0.0, cdist(track.smooth_feat.reshape(1,-1), det_features, metric))
    track_features = np.asarray([track.smooth_feat for track in tracks], dtype=np.float32)
    if metric == "cosine":
        track_norm = np.linalg.norm(track_features, axis=1, keepdims=True)
        det_norm = np.linalg.norm(det_features, axis=1, keepdims=True).T
        cost_matrix = 1 - track_features @ det_features.T / np.maximum(track_norm * det_norm, np.finfo(float).eps)
    else:
        from scipy.spatial.distance import cdist

        cost_matrix = cdist(track_features, det_features, metric)
    cost_matrix = np.maximum(0.0, cost_matrix)  # Normalized features
    return cost_matrix