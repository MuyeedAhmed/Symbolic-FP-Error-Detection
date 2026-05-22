    def grad_hess(x):
        return grad(x), lambda x: A.T.dot(A.dot(x))