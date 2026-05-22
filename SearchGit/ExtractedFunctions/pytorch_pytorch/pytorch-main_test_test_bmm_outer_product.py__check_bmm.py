    def _check_bmm(self, a, b, **kwargs):
        self.assertEqual(torch.bmm(a, b), a @ b, **kwargs)