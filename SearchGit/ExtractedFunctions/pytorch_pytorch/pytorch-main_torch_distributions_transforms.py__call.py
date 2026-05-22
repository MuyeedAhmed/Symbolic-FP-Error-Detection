    def _call(self, x):
        x = LowerCholeskyTransform()(x)
        return x @ x.mT