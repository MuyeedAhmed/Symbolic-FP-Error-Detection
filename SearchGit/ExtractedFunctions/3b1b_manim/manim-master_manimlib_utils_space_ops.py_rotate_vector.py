def rotate_vector(
    vector: Vect3,
    angle: float,
    axis: Vect3 = OUT
) -> Vect3:
    rot = Rotation.from_rotvec(angle * normalize(axis))
    return np.dot(vector, rot.as_matrix().T)