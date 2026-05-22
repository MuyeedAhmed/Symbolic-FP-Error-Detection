    def test_dtensor_mm(self):
        """
        Test mm with DTensor with 2D mesh.
        We need to add the test here since we only test 1D mesh in test_dtensor_ops.py.
        Also, we added tests for the corner case where one of the 2D dimension is 1.

        # TODO: we need to test more DTensor ops with 2D mesh, especially when 1 of the
        mesh dimension of the 2D mesh is 1.
        """
        mesh_0 = init_device_mesh(self.device_type, (self.world_size // 2, 2))
        mesh_1 = init_device_mesh(self.device_type, (self.world_size, 1))
        mesh_2 = init_device_mesh(self.device_type, (1, self.world_size))

        for mesh in [mesh_0, mesh_1, mesh_2]:
            lhs = torch.randn(256, 128)
            rhs = torch.randn(128, 256)
            mm_result = lhs @ rhs

            lhs_dtensor = distribute_tensor(lhs, mesh, [Shard(dim=0), Replicate()])
            rhs_dtensor = distribute_tensor(rhs, mesh, [Replicate(), Shard(dim=1)])
            dtensor_result = lhs_dtensor @ rhs_dtensor
            self.assertEqual(
                dtensor_result.full_tensor(), mm_result, atol=1.5e-5, rtol=1e-6
            )