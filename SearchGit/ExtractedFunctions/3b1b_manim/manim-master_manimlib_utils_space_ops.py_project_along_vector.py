def project_along_vector(point: Vect3, vector: Vect3) -> Vect3:
    matrix = np.identity(3) - np.outer(vector, vector)
    return np.dot(point, matrix.T)