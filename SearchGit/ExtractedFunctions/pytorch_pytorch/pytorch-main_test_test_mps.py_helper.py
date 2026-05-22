        def helper(n, atol=1e-5, rtol=1e-6):
            # Generate well-conditioned invertible matrix by adding scaled identity
            # This ensures the matrix is not singular
            cpu_input = torch.randn(n, n, device='cpu') + torch.eye(n, device='cpu') * 10
            mps_input = cpu_input.to('mps')

            cpu_result = torch.linalg.inv(cpu_input)
            mps_result = torch.linalg.inv(mps_input)
            self.assertEqual(cpu_result, mps_result, atol=atol, rtol=rtol)