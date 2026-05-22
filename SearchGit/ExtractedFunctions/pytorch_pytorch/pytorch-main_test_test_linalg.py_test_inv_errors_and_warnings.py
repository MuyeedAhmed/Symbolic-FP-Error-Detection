    def test_inv_errors_and_warnings(self, device, dtype):
        # inv expects batches of square matrices as input
        a = torch.randn(2, 3, 4, 3, dtype=dtype, device=device)
        with self.assertRaisesRegex(RuntimeError, "must be batches of square matrices"):
            torch.linalg.inv(a)

        # inv requires the input to be at least 2 dimensional tensor
        a = torch.randn(2, device=device, dtype=dtype)
        with self.assertRaisesRegex(RuntimeError, "must have at least 2 dimensions"):
            torch.linalg.inv(a)

        # if input is not invertible, RuntimeError is raised mentioning the first non-invertible batch
        def run_test_singular_input(batch_dim, n):
            a = torch.eye(3, 3, dtype=dtype, device=device).reshape((1, 3, 3)).repeat(batch_dim, 1, 1)
            a[n, -1, -1] = 0
            with self.assertRaisesRegex(torch.linalg.LinAlgError, rf"\(Batch element {n}\): The diagonal element 3 is zero"):
                torch.linalg.inv(a)

        for params in [(1, 0), (2, 0), (2, 1), (4, 0), (4, 2), (10, 2)]:
            run_test_singular_input(*params)

        # dtypes should match
        a = torch.eye(2, dtype=dtype, device=device)
        out = torch.empty(0, dtype=torch.int, device=device)
        with self.assertRaisesRegex(RuntimeError, "but got int instead"):
            torch.linalg.inv(a, out=out)

        # device should match
        if torch.cuda.is_available():
            wrong_device = 'cpu' if self.device_type != 'cpu' else 'cuda'
            out = torch.empty(0, device=wrong_device, dtype=dtype)
            with self.assertRaisesRegex(RuntimeError, "tensors to be on the same device"):
                torch.linalg.inv(a, out=out)

        # if out tensor with wrong shape is passed a warning is given
        with warnings.catch_warnings(record=True) as w:
            a = torch.eye(2, dtype=dtype, device=device)
            out = torch.empty(1, dtype=dtype, device=device)
            # Trigger warning
            torch.linalg.inv(a, out=out)
            # Check warning occurs
            self.assertEqual(len(w), 1)
            self.assertTrue("An output with one or more elements was resized" in str(w[-1].message))

        # if out tensor in batched column major format but with wrong a warning is given
        with warnings.catch_warnings(record=True) as w:
            a = torch.eye(2, dtype=dtype, device=device)
            out = torch.empty(3, 3, dtype=dtype, device=device)
            out = out.mT.clone(memory_format=torch.contiguous_format)
            out = out.mT
            self.assertTrue(out.mT.is_contiguous())
            # Trigger warning
            torch.linalg.inv(a, out=out)
            # Check warning occurs
            self.assertEqual(len(w), 1)
            self.assertTrue("An output with one or more elements was resized" in str(w[-1].message))