        def run_test_singular_input(batch_dim, n):
            a = torch.eye(3, 3, dtype=dtype, device=device).reshape((1, 3, 3)).repeat(batch_dim, 1, 1)
            a[n, -1, -1] = 0
            with self.assertRaisesRegex(torch.linalg.LinAlgError, rf"\(Batch element {n}\): The diagonal element 3 is zero"):
                torch.linalg.inv(a)