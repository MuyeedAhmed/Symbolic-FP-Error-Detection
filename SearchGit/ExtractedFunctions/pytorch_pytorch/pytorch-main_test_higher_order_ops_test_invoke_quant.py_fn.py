        def fn(x, y, z):
            return (
                invoke_quant_tracer(
                    gn, x, y, scheme="nf4", quant_options=invoke_quant_tracer
                )
                @ z
            )