    def test_addmm_badmm_scalar_tnesor_input(self, device, dtype):
        input = torch.tensor(1).to(device=device, dtype=dtype)

        # test addmm
        mat1 = torch.randn(10, 25, device=device).to(dtype)
        mat2 = torch.randn(25, 10, device=device).to(dtype)
        result = torch.addmm(input, mat1, mat2)

        ref = mat1.cpu().numpy() @ mat2.cpu().numpy() + 1
        self.assertEqual(result, ref)

        # test baddbmm
        mat1 = torch.randn(3, 10, 25, device=device).to(dtype)
        mat2 = torch.randn(3, 25, 10, device=device).to(dtype)
        result = torch.baddbmm(input, mat1, mat2)

        ref = mat1.cpu().numpy() @ mat2.cpu().numpy() + 1
        self.assertEqual(result, ref)