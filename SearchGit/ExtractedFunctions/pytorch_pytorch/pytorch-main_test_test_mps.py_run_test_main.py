        def run_test_main(A, hermitian):
            # Testing against definition for pseudo-inverses
            A_pinv = torch.linalg.pinv(A, hermitian=hermitian)
            np_A = A.cpu().numpy()
            np_A_pinv = A_pinv.cpu().numpy()
            if A.numel() > 0:
                self.assertEqual(A, np_A @ np_A_pinv @ np_A, atol=precision, rtol=precision)
                self.assertEqual(A_pinv, np_A_pinv @ np_A @ np_A_pinv, atol=precision, rtol=precision)
                self.assertEqual(np_A @ np_A_pinv, (np_A @ np_A_pinv).conj().swapaxes(-2, -1), atol=precision, rtol=precision)
                self.assertEqual(np_A_pinv @ np_A, (np_A_pinv @ np_A).conj().swapaxes(-2, -1), atol=precision, rtol=precision)
            else:
                self.assertEqual(A.shape, A_pinv.shape[:-2] + (A_pinv.shape[-1], A_pinv.shape[-2]))

            # Check out= variant
            out = torch.empty_like(A_pinv)
            ans = torch.linalg.pinv(A, hermitian=hermitian, out=out)
            self.assertEqual(ans, out)
            self.assertEqual(ans, A_pinv)