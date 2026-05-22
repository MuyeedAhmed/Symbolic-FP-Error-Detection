    def test_linalg_ops(self):
        device_mesh = DeviceMesh(self.device_type, list(range(self.world_size)))

        # Batched positive-definite matrix for ops that require it (cholesky, inv)
        A = torch.randn(8, 4, 4, device=self.device_type)
        spd = A @ A.mT + 4 * torch.eye(4, device=self.device_type)
        dt_spd = distribute_tensor(spd, device_mesh, [Shard(0)])

        # cholesky
        expected = torch.linalg.cholesky(spd)
        result = torch.linalg.cholesky(dt_spd)
        self.assertEqual(result.full_tensor(), expected)
        self.assertTrue(result.placements[0].is_shard(0))

        # inv
        expected = torch.linalg.inv(spd)
        result = torch.linalg.inv(dt_spd)
        self.assertEqual(result.full_tensor(), expected)
        self.assertTrue(result.placements[0].is_shard(0))

        # lu_factor
        B = torch.randn(8, 4, 4, device=self.device_type)
        dt_B = distribute_tensor(B, device_mesh, [Shard(0)])
        exp_LU, exp_piv = torch.linalg.lu_factor(B)
        res_LU, res_piv = torch.linalg.lu_factor(dt_B)
        self.assertEqual(res_LU.full_tensor(), exp_LU)
        self.assertTrue(res_LU.placements[0].is_shard(0))

        # eig
        expected_vals, expected_vecs = torch.linalg.eig(B)
        result_vals, result_vecs = torch.linalg.eig(dt_B)
        self.assertEqual(result_vals.full_tensor(), expected_vals)
        self.assertTrue(result_vals.placements[0].is_shard(0))

        # solve
        rhs = torch.randn(8, 4, 2, device=self.device_type)
        dt_rhs = distribute_tensor(rhs, device_mesh, [Shard(0)])
        expected = torch.linalg.solve(spd, rhs)
        result = torch.linalg.solve(dt_spd, dt_rhs)
        self.assertEqual(result.full_tensor(), expected)
        self.assertTrue(result.placements[0].is_shard(0))