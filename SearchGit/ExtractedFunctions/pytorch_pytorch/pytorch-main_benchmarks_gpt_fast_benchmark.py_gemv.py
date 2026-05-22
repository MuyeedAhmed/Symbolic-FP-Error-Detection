            def gemv(W, x):
                return W.to(x.dtype) @ x