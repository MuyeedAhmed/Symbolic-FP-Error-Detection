    def test_profile_loading_and_lookup(self):
        """Load a trace with collectives, aten ops, and a custom op; verify lookups."""
        with torch.library._scoped_library("test_pge", "DEF") as lib:
            lib.define("my_op(Tensor x, Tensor w) -> Tensor")
            lib.impl("my_op", lambda x, w: x @ w, "CPU")
            lib.impl("my_op", lambda x, w: x @ w, "Meta")

            trace = _make_pge_trace(
                collectives=[
                    {
                        "name": "allreduce",
                        "dur": 100.0,
                        "nelems": 1000,
                        "dtype": "Float",
                        "ranks": "[0,2,4,6]",
                        "group_size": 4,
                    },
                    {
                        "name": "allreduce",
                        "dur": 800.0,
                        "nelems": 8000,
                        "dtype": "Float",
                        "ranks": "[0,2,4,6]",
                        "group_size": 4,
                    },
                ],
                matmuls=[
                    {
                        "shapes": [[128, 256], [256, 512]],
                        "dur": 50.0,
                        "dtypes": ["float", "float"],
                    }
                ],
                pg_config={"0": {"ranks": [0, 2, 4, 6]}},
            )
            # Inject a custom op event
            eid = 9000
            trace["traceEvents"].extend(
                [
                    {
                        "cat": "cpu_op",
                        "name": "test_pge::my_op",
                        "dur": 0,
                        "args": {
                            "External id": eid,
                            "Input Dims": [[32, 64], [64, 128]],
                            "Input Strides": [[64, 1], [128, 1]],
                            "Input type": ["c10::BFloat16", "c10::BFloat16"],
                        },
                    },
                    {
                        "cat": "kernel",
                        "dur": 42.0,
                        "name": "custom_kernel",
                        "args": {"External id": eid},
                    },
                ]
            )
            profile = _load_pge_profile(trace)

            # Collectives
            self.assertAlmostEqual(
                profile.lookup_collective("all_reduce", (0, 2, 4, 6), 1000, "Float")[0],
                0.1,
                places=4,
            )
            self.assertGreater(
                profile.lookup_collective("all_reduce", (0, 2, 4, 6), 4000, "Float")[0],
                0.1,
            )
            self.assertIsNotNone(
                profile.lookup_collective("all_reduce", (1, 3, 5, 7), 1000, "Float")
            )
            self.assertIsNone(
                profile.lookup_collective("all_to_all", (0, 2, 4, 6), 1000, "Float")
            )

            # Aten op
            self.assertAlmostEqual(
                profile.lookup_op(
                    "aten::mm",
                    ((128, 256), (256, 512)),
                    ((256, 1), (512, 1)),
                    torch.float32,
                ),
                0.05,
                places=4,
            )
            self.assertIsNone(
                profile.lookup_op(
                    "aten::mm", ((999,), (999,)), ((1,), (1,)), torch.float32
                )
            )

            # Custom op
            self.assertAlmostEqual(
                profile.lookup_op(
                    "test_pge::my_op",
                    ((32, 64), (64, 128)),
                    ((64, 1), (128, 1)),
                    torch.bfloat16,
                ),
                0.042,
                places=4,
            )

            # End-to-end: FX graph with custom op -> estimator match
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".json", delete=False
            ) as f:
                json.dump(trace, f)
                f.flush()
                trace_path = f.name
            try:
                estimator = ProfileGuidedEstimator(trace_path)
                with FakeTensorMode():
                    x = torch.randn(32, 64, dtype=torch.bfloat16, device="cpu")
                    w = torch.randn(64, 128, dtype=torch.bfloat16, device="cpu")
                    gm = make_fx(torch.ops.test_pge.my_op)(x, w)
                self.assertTrue(
                    any(estimator(n) is not None for n in gm.graph.nodes),
                    "Estimator should match custom op",
                )
            finally:
                os.unlink(trace_path)