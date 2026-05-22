    def test_1_sized_with_0_strided(self, device, dtype):
        a = make_tensor((8, 1, 64), dtype=dtype, device=device)
        a_strided = torch.as_strided(a, size=[8, 1, 64], stride=[64, 0, 1])
        b = make_tensor((8, 64, 512), dtype=dtype, device=device)
        b_strided = torch.as_strided(b, size=[8, 64, 512], stride=[64, 1, 512])
        res = torch.bmm(a_strided, b_strided)
        expect = torch.from_numpy(a_strided.cpu().numpy() @ b_strided.cpu().numpy()).to(
            device=device, dtype=dtype
        )
        self.assertEqual(expect, res)