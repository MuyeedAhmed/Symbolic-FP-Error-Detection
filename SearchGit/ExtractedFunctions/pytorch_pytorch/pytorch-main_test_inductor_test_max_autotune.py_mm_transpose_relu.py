        def mm_transpose_relu(a, b):
            return (a @ b.transpose(0, 1)).relu()