    def test_mm_with_mH_args(self, dtype, backend):
        # Testing mm with mH-transformed arguments.
        # The root cause of:
        # https://github.com/pytorch/pytorch/issues/174382
        val = 3 + 4j
        x = torch.zeros(2, 3, dtype=dtype, device="cuda")
        x.diagonal().fill_(val)

        ref_corrcoef = torch.empty(2, 2, dtype=dtype, device="cuda")
        ref_corrcoef.fill_(-0.5)
        ref_corrcoef.diagonal().fill_(1.0)

        with blas_library_context(backend):
            for a in (x, x.mH):
                norm_squared = (a @ a.mH).sum().item()
                self.assertEqual(norm_squared, 50 + 0j)

            self.assertEqual(torch.corrcoef(x), ref_corrcoef)