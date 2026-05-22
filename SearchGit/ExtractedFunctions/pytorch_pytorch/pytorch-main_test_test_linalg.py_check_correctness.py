        def check_correctness(a, b):
            sol = torch.linalg.lstsq(a, b).solution
            sol2 = a.pinverse() @ b
            self.assertEqual(sol, sol2, rtol=1e-5, atol=1e-5)