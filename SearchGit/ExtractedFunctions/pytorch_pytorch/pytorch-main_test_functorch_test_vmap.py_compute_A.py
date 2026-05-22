        def compute_A(out):
            L, Q = out
            n = Q.shape[-1]
            diag_L = L.new_zeros(*L.shape[:-1], n, n)
            diag_L.diagonal(offset=0, dim1=-2, dim2=-1).copy_(L)
            Qh = Q.transpose(-2, -1).conj()
            return Q @ diag_L @ Qh