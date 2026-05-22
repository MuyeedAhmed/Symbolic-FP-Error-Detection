    def test_matmul_device_mismatch(self):
        cpu = torch.rand((10, 10))
        cuda = cpu.cuda()
        with self.assertRaisesRegex(
            RuntimeError, "Expected all tensors to be on the same device"
        ):
            cpu @ cuda
        with self.assertRaisesRegex(
            RuntimeError, "Expected all tensors to be on the same device"
        ):
            cuda @ cpu

        for s, m1, m2 in product((cpu, cuda), repeat=3):
            if s.device == m1.device == m2.device:
                torch.addmm(s, m1, m2)
            else:
                with self.assertRaisesRegex(
                    RuntimeError, "Expected all tensors to be on the same device"
                ):
                    torch.addmm(s, m1, m2)