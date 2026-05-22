    def test_sampled_addmm_autograd(self, device, dtype):
        from torch.testing._internal.common_methods_invocations import sample_inputs_sparse_sampled_addmm

        samples = list(sample_inputs_sparse_sampled_addmm(None, device, dtype, requires_grad=True))

        for sample, dense_covector in zip(samples, [True, False]):
            c = sample.input
            a = sample.args[0]
            b = sample.args[1]

            # Compute sparse result
            output = torch.sparse.sampled_addmm(c, a, b, **sample.kwargs)
            covector = torch.randn_like(output).to_dense() if dense_covector else torch.randn_like(output)
            output.backward(covector)

            # Compute dense result and compare with sparse result
            c1, a1, b1 = (x.detach().to_dense().requires_grad_(True) for x in [c, a, b])
            dense_output = sample.kwargs['alpha'] * (a1 @ b1) * torch.ones_like(c).to_dense() + sample.kwargs['beta'] * c1
            self.assertEqual(output, dense_output)
            dense_covector = covector.to_dense()
            dense_output.backward(dense_covector)
            self.assertEqual(c.grad, c1.grad)
            self.assertEqual(a.grad, a1.grad)
            self.assertEqual(b.grad, b1.grad)