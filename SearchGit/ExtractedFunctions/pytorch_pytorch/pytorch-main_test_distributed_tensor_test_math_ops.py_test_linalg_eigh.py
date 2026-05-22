    def test_linalg_eigh(self):
        A = torch.randn(2, 2, dtype=torch.float64)
        mesh = self.build_device_mesh()
        dtensor_A = distribute_tensor(A, device_mesh=mesh, placements=[Replicate()])
        dtensor_A = dtensor_A + dtensor_A.mT
        dtensor_L, dtensor_Q = torch.linalg.eigh(dtensor_A)

        # TODO: we need to convert A, L, Q to local because we don't have a
        # sharding strategy registered for aten.dist.default yet.
        local_A, local_L, local_Q = (
            dtensor_A.to_local(),
            dtensor_L.to_local(),
            dtensor_Q.to_local(),
        )
        distance = torch.dist(local_Q @ torch.diag(local_L) @ local_Q.mT, local_A)
        self.assertEqual(distance.item(), 0.0)