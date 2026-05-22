            def gather_gemv(W, score_idxs, x):
                return W[score_idxs].to(x.dtype) @ x