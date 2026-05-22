        def ref_sparse_mm(a, b):
            return a.to_dense() @ b.to_dense()