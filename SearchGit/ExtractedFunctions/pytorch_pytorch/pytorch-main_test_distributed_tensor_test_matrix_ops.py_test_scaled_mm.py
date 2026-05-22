    def test_scaled_mm(self):
        device_mesh = self.build_device_mesh()
        shrd0 = Shard(0)
        shrd1 = Shard(1)
        repl = Replicate()
        part = Partial()

        ws = self.world_size
        # _scaled_mm requires all dimensions to be multiples of 16. Since we'll
        # shard along n and k, we need to ensure this stays true on each rank.
        m, n, k = 16, 32 * ws, 16 * ws

        t1 = torch.randn(m, k, device=self.device_type, dtype=torch.bfloat16)
        t2 = torch.randn(n, k, device=self.device_type, dtype=torch.bfloat16)

        for (
            output_spec,
            t1_spec,
            t2_spec,
            scale1_shape,
            scale2_shape,
            scale1_spec,
            scale2_spec,
        ) in [
            # Tensor-wise scaling
            # Replicated, zero-dim scale
            (repl, repl, repl, (), (), repl, repl),
            # Column-parallel, two-dim scale
            (shrd1, repl, shrd0, (1, 1), (1, 1), repl, repl),
            # Row-parallel, one-dim scale
            (part, shrd1, shrd1, (1,), (1,), repl, repl),
            # Row-wise scaling
            # Replicated
            (repl, repl, repl, (m, 1), (n, 1), repl, repl),
            # Column-parallel
            (shrd1, repl, shrd0, (m, 1), (n, 1), repl, shrd0),
            # Row-parallel (which actually ends up doing sub-row-wise scaling)
            (part, shrd1, shrd1, (m, ws), (n, ws), shrd1, shrd1),
        ]:
            full_ref_res = t1 @ t2.t()

            t1_fp8, scale1 = scale_for_fp8(t1, scale1_shape)
            t2_fp8, scale2 = scale_for_fp8(t2, scale2_shape)

            dist_t1_fp8 = distribute_tensor(t1_fp8, device_mesh, [t1_spec])
            dist_t2_fp8 = distribute_tensor(t2_fp8, device_mesh, [t2_spec])
            dist_scale1 = distribute_tensor(scale1, device_mesh, [scale1_spec])
            dist_scale2 = distribute_tensor(scale2, device_mesh, [scale2_spec])

            with CommDebugMode() as comm_mode:
                dist_res = cast(
                    DTensor,
                    torch._scaled_mm(
                        dist_t1_fp8,
                        dist_t2_fp8.t(),
                        scale_a=dist_scale1,
                        scale_b=dist_scale2.t(),
                        out_dtype=torch.bfloat16,
                    ),
                )

            self.assertEqual(dist_res.placements[0], output_spec)

            full_dist_res = dist_res.full_tensor()
            # Fp8 matmuls are quite inaccurate, we need high tolerances
            self.assertEqual(full_dist_res, full_ref_res, atol=1.5, rtol=7e-2)

            self.assertEqual(comm_mode.get_total_counts(), 0)