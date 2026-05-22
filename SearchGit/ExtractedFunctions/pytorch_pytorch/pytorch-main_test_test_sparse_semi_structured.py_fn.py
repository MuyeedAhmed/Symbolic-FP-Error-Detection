        def fn(x):
            y = SparseSemiStructuredTensorCUSPARSELT.prune_dense_static_sort(x)
            y = y.t()
            return x @ y