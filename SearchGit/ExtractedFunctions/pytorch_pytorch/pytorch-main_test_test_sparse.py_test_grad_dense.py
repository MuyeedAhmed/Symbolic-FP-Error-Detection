            def test_grad_dense(a_s, b_s, g_s):
                a = a_s.to_dense().detach()
                b = b_s.to_dense().detach()
                g = g_s.to_dense().detach()

                a.requires_grad_(True)
                b.requires_grad_(True)
                c = a @ b
                c.backward(g)
                return a.grad.sparse_mask(a_s.coalesce()), b.grad.sparse_mask(b_s.coalesce())