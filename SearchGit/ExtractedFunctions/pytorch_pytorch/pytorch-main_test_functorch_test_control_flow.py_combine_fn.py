        def combine_fn(carry, xs):
            # carry: 8, 5
            # xs: 5, 8
            # new_carry: 8, 5
            new_carry = carry + (xs * scale).sum()
            output = xs @ carry
            return new_carry, output