        def fn(a, b):
            # Chain of 4 matmuls: mm0 -> mm1 -> mm2 -> mm3
            mm0 = a @ b
            mm1 = mm0 @ b
            mm2 = mm1 @ b
            mm3 = mm2 @ b
            return mm3