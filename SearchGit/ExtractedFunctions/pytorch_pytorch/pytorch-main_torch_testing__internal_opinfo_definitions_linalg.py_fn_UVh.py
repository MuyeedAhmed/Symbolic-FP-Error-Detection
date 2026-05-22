    def fn_UVh(usv):
        U, S, Vh = uniformize(usv)
        return U @ Vh, S