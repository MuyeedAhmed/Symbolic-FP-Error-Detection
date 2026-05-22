        def fn(vec):
            # Expand inside the function to make sure the input to
            # gradcheck does not have overlapping memory
            vec = vec.expand(2)
            return (mat @ vec).sum()