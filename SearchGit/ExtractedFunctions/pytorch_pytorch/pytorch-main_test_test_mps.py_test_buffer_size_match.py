    def test_buffer_size_match(self):
        # this test shouldn't cause any crash
        size = 16
        cpu_A = torch.rand(size, device='cpu')
        cpu_F = torch.rand(size, size, size, device='cpu')

        mps_A = cpu_A.to('mps')
        mps_F = cpu_F.to('mps')
        self.assertEqual(cpu_A @ cpu_F, mps_A @ mps_F)