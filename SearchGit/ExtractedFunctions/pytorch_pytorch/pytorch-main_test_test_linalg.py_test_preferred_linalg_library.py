    def test_preferred_linalg_library(self):
        # The main purpose of this test is to make sure these "backend" calls work normally without raising exceptions.
        x = torch.randint(2, 5, (2, 4, 4), device='cuda', dtype=torch.double)

        torch.backends.cuda.preferred_linalg_library('cusolver')
        out1 = torch.linalg.inv(x)

        torch.backends.cuda.preferred_linalg_library('default')
        # Although linalg preferred flags doesn't affect CPU currently,
        # we set this to make sure the flag can switch back to default normally.
        out_ref = torch.linalg.inv(x.cpu())

        self.assertEqual(out_ref, out1.cpu())

        if torch.cuda.has_magma:
            torch.backends.cuda.preferred_linalg_library('magma')
            out2 = torch.linalg.inv(x)
            self.assertEqual(out1, out2)