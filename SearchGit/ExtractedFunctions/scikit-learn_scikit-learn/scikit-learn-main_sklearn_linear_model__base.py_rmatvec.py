                def rmatvec(b):
                    return X.T.dot(b) - X_offset * b.sum()