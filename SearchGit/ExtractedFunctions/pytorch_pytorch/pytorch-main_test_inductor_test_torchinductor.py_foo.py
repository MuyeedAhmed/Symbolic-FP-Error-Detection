        def foo(inp, weight):
            matmul_output = inp @ weight
            final_output = layer_norm(matmul_output)
            return final_output