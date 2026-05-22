    def fn(usv):
        U, S, V = usv
        return U @ V.mH, S