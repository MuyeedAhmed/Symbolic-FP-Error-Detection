        def tril_cholesky_to_tril_corr(x):
            x = vec_to_tril_matrix(x, -1)
            diag = (1 - (x * x).sum(-1)).sqrt().diag_embed()
            x = x + diag
            return tril_matrix_to_vec(x @ x.T, -1)