            def combine_fn(carry, x):
                new_carry = {
                    "param": carry["param"] @ x + carry["bias"],
                    "bias": carry["bias"].sin(),
                }
                return new_carry, (
                    pytree.tree_map(lambda x: x.clone(), new_carry),
                    {"dummy": x.sin()},
                )