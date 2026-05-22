        def mul_svd_factors(U, S, Vh):
            return (U * S.to(dtype).unsqueeze(-2)) @ Vh