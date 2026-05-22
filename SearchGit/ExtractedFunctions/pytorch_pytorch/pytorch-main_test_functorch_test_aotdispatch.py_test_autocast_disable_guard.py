    def test_autocast_disable_guard(self):
        with torch._C._DisableAutocast():
            x = torch.rand([4, 4]).cuda()
            y = x @ x
            self.assertEqual(y.dtype, torch.float32)