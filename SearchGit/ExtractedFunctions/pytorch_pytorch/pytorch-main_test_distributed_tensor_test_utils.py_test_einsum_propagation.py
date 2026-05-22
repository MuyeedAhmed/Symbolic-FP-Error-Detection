    def test_einsum_propagation(self):
        with LocalTensorMode(ranks=self.world_size):
            mesh = init_device_mesh("cpu", (4, 4))
            input_tensor = torch.arange(16 * 16).float().view(16, 16)
            A = distribute_tensor(input_tensor, mesh, [Shard(1), Shard(1)])
            B1 = distribute_tensor(input_tensor, mesh, [Shard(0), Shard(0)])
            B2 = distribute_tensor(
                input_tensor,
                mesh,
                [_StridedShard(0, split_factor=mesh.size(1)), Shard(0)],
            )
            with CommDebugMode() as comm_mode:
                # res1 will be (Partial, Partial), no redistribution needed
                res1 = A @ B1
            self.assertEqual(
                comm_mode.get_comm_counts()[c10d_functional.all_gather_into_tensor], 0
            )

            with CommDebugMode() as comm_mode:
                # `A @ B2` will trigger redistribution on both inputs as below:
                # A: S(1)[0]S(1)[1]->S(1)R->RR
                # B2: S(0)[1]S(0)[0]->RS(0)->RR
                res2 = A @ B2
            self.assertEqual(
                comm_mode.get_comm_counts()[c10d_functional.all_gather_into_tensor], 4
            )
            if not isinstance(res1, DTensor):
                raise AssertionError(f"Expected res1 to be DTensor, got {type(res1)}")
            if not isinstance(res2, DTensor):
                raise AssertionError(f"Expected res2 to be DTensor, got {type(res2)}")
            self.assertEqual(res1.full_tensor(), res2.full_tensor())