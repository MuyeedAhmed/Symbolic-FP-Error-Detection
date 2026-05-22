        def model(a, y):
            b = a + 1
            b2 = b - 1
            c = y @ y  # mm boundary
            d = b * 2
            d2 = d / 2
            return b2, c, d2