    def test_eig_cuda_complex_eigenvectors(self, device, dtype):
        """Test CUDA eigenvector decoding with known ground truth, including batching."""

        # Test 1: Rotation matrix (complex eigenvalues - conjugate pairs)
        theta = math.pi / 4
        A_complex = torch.tensor([
            [math.cos(theta), -math.sin(theta)],
            [math.sin(theta), math.cos(theta)]
        ], dtype=dtype, device=device)

        vals_complex, vecs_complex = torch.linalg.eig(A_complex)

        # Verify eigenvalues are e^(±iθ) for rotation by θ
        # For θ = π/4, eigenvalues are e^(±iπ/4) - a conjugate pair
        expected_eigenvalue = complex(math.cos(theta), math.sin(theta))
        expected_val = torch.tensor(
            expected_eigenvalue, dtype=vals_complex.dtype, device=device
        )
        expected_val_conj = torch.tensor(
            expected_eigenvalue.conjugate(), dtype=vals_complex.dtype, device=device
        )
        # Check both eigenvalues are present and form a conjugate pair
        match_0_pos = torch.allclose(vals_complex[0], expected_val, atol=1e-5, rtol=1e-5)
        match_0_neg = torch.allclose(vals_complex[0], expected_val_conj, atol=1e-5, rtol=1e-5)
        match_1_pos = torch.allclose(vals_complex[1], expected_val, atol=1e-5, rtol=1e-5)
        match_1_neg = torch.allclose(vals_complex[1], expected_val_conj, atol=1e-5, rtol=1e-5)
        # Valid if (vals[0]=λ AND vals[1]=λ*) OR (vals[0]=λ* AND vals[1]=λ)
        self.assertTrue(
            (match_0_pos and match_1_neg) or (match_0_neg and match_1_pos),
            f"Expected conjugate pair {{λ, λ*}}, got {vals_complex[0]}, {vals_complex[1]}"
        )

        # Verify output is complex type
        self.assertTrue(vals_complex.dtype in [torch.complex64, torch.complex128])
        self.assertTrue(vecs_complex.dtype in [torch.complex64, torch.complex128])

        # Verify Av = λv for all eigenpairs (vectorized)
        lhs = A_complex.to(vecs_complex.dtype) @ vecs_complex
        rhs = vals_complex.unsqueeze(-2) * vecs_complex
        self.assertEqual(lhs, rhs, atol=1e-5, rtol=1e-5)

        # Test 2: Diagonal matrix (all real eigenvalues)
        A_real = torch.diag(torch.tensor([1.0, 2.0, 3.0], dtype=dtype, device=device))

        vals_real, vecs_real = torch.linalg.eig(A_real)

        # Output is still complex type, but imaginary parts should be ~zero
        self.assertTrue(torch.allclose(vals_real.imag, torch.zeros_like(vals_real.imag), atol=1e-6))
        # Real parts should match diagonal values
        self.assertTrue(torch.allclose(
            torch.sort(vals_real.real)[0],
            torch.tensor([1., 2., 3.], dtype=dtype, device=device),
            atol=1e-6, rtol=1e-6
        ))

        # Verify Av = λv for all eigenpairs (vectorized)
        lhs = A_real.to(vecs_real.dtype) @ vecs_real
        rhs = vals_real.unsqueeze(-2) * vecs_real
        self.assertEqual(lhs, rhs, atol=1e-5, rtol=1e-5)

        # Test 3: Batched - mix of real and complex eigenvalues
        A_batch = torch.stack([
            # Rotation (complex eigenvalues)
            torch.tensor([
                [math.cos(math.pi / 6), -math.sin(math.pi / 6)],
                [math.sin(math.pi / 6), math.cos(math.pi / 6)]
            ], dtype=dtype, device=device),
            # Diagonal (real eigenvalues)
            torch.diag(torch.tensor([4.0, 5.0], dtype=dtype, device=device)),
            # Another rotation (complex eigenvalues)
            torch.tensor([
                [math.cos(math.pi / 3), -math.sin(math.pi / 3)],
                [math.sin(math.pi / 3), math.cos(math.pi / 3)]
            ], dtype=dtype, device=device),
        ])

        vals_batch, vecs_batch = torch.linalg.eig(A_batch)

        # Verify Av = λv for all matrices in batch
        lhs = A_batch.to(vecs_batch.dtype) @ vecs_batch
        rhs = vals_batch.unsqueeze(-2) * vecs_batch
        self.assertEqual(lhs, rhs, atol=1e-5, rtol=1e-5)