    def matvec(b):
        return X.dot(b) - sample_weight_sqrt * b.dot(X_offset)