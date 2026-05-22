    def combine_fn(carry, x):
        result = carry @ x + x
        return result, carry.clone()