def angle_between_vectors(v1: VectN, v2: VectN) -> float:
    """
    Returns the angle between two 3D vectors.
    This angle will always be btw 0 and pi
    """
    n1 = get_norm(v1)
    n2 = get_norm(v2)
    if n1 == 0 or n2 == 0:
        return 0
    cos_angle = np.dot(v1, v2) / np.float64(n1 * n2)
    return math.acos(clip(cos_angle, -1, 1))