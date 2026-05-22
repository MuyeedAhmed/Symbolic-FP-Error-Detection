            def reference_model(a, b, bias):
                matmul_result = a @ b
                biased = matmul_result + bias
                activated = torch.relu(biased)
                scaled = activated * 2.0
                return scaled