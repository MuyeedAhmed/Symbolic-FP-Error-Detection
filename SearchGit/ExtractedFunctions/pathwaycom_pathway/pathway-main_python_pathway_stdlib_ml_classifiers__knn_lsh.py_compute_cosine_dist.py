def compute_cosine_dist(data_table: np.ndarray, query_point: np.ndarray) -> np.ndarray:
    return 1 - np.dot(data_table, query_point) / (
        np.linalg.norm(data_table, axis=1) * np.linalg.norm(query_point)
    )