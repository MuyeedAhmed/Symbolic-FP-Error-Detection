def scipy_matmul(mat1, mat2):
    if isspmatrix(mat1) and isspmatrix(mat2):
        return mat1.dot(mat2).tocoo()
    return mat1.dot(mat2)