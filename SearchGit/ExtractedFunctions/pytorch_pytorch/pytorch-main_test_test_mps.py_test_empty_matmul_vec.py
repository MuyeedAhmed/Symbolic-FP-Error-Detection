    def test_empty_matmul_vec(self):
        tensor_1 = torch.rand((0, 100), device="mps")
        tensor_2 = torch.rand((100, ), device="mps")
        self.assertEqual((tensor_1 @ tensor_2).cpu(), tensor_1.cpu() @ tensor_2.cpu())