    def test_linalg_solve_partial(self):
        device_mesh = DeviceMesh(self.device_type, list(range(self.world_size)))
        A = torch.randn(4, 4, device=self.device_type, dtype=torch.float64)
        A = A @ A.mT + 4 * torch.eye(4, device=self.device_type, dtype=torch.float64)
        B = torch.randn(4, 2, device=self.device_type, dtype=torch.float64)

        dt_A = distribute_tensor(A, device_mesh, [Replicate()])
        dt_B = distribute_tensor(B, device_mesh, [Partial()])
        expected = torch.linalg.solve(A, B)
        result = torch.linalg.solve(dt_A, dt_B)
        self.assertEqual(result.full_tensor(), expected)