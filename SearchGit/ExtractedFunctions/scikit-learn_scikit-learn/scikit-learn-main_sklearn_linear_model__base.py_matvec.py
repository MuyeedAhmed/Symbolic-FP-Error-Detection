                def matvec(b):
                    return X.dot(b) - b.dot(X_offset)