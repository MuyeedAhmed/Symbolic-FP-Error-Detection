    def test_bmm(self, device, dtype):
        batch_sizes = [1, 10]
        M, N, O = 23, 15, 12
        numpy_dtype = dtype if dtype != torch.bfloat16 else torch.float32

        def invert_perm(p):
            d = {x: i for i, x in enumerate(p)}
            return (d[0], d[1], d[2])

        def generate_inputs(num_batches):
            # transposed tensors
            for perm1, perm2 in itertools.product(
                itertools.permutations((0, 1, 2)), repeat=2
            ):
                b1 = make_tensor(
                    (num_batches, M, N), dtype=dtype, device=device, low=-0.1, high=0.1
                )
                b2 = make_tensor(
                    (num_batches, N, O), dtype=dtype, device=device, low=-0.1, high=0.1
                )
                b1 = b1.permute(perm1).contiguous().permute(invert_perm(perm1))
                b2 = b2.permute(perm2).contiguous().permute(invert_perm(perm2))
                yield b1, b2
            # broadcasting tensors
            for b1, b2, b3, b4, b5, b6 in itertools.product((True, False), repeat=6):
                shape1 = (num_batches if b1 else 1, M if b2 else 1, N if b3 else 1)
                shape2 = (num_batches if b4 else 1, N if b5 else 1, O if b6 else 1)
                b1 = make_tensor(
                    shape1, dtype=dtype, device=device, low=-0.1, high=0.1
                ).expand(num_batches, M, N)
                b2 = make_tensor(
                    shape2, dtype=dtype, device=device, low=-0.1, high=0.1
                ).expand(num_batches, N, O)
                yield b1, b2
            # zero-sized tensors
            for z1, z2, z3, z4 in itertools.product((True, False), repeat=4):
                shape1 = (num_batches if z1 else 0, M if z2 else 0, N if z3 else 0)
                shape2 = (num_batches if z1 else 0, N if z3 else 0, O if z4 else 0)
                b1 = torch.randn(shape1, dtype=dtype, device=device)
                b2 = torch.randn(shape2, dtype=dtype, device=device)
                yield b1, b2

        for num_batches in batch_sizes:
            for (b1, b2), perm3 in itertools.product(
                generate_inputs(num_batches), itertools.permutations((0, 1, 2))
            ):
                res1 = torch.bmm(b1, b2)
                res2 = (
                    torch.full(
                        (num_batches, M, O), math.nan, dtype=dtype, device=device
                    )
                    .permute(perm3)
                    .contiguous()
                    .permute(invert_perm(perm3))
                )
                torch.bmm(b1, b2, out=res2)
                expect = torch.from_numpy(
                    b1.to(numpy_dtype).cpu().numpy() @ b2.to(numpy_dtype).cpu().numpy()
                ).to(device=device, dtype=dtype)
                self.assertEqual(expect, res1)
                self.assertEqual(expect, res2)

                if self.device_type == "cuda":
                    # check that mixed arguments are rejected
                    self.assertRaises(RuntimeError, lambda: torch.bmm(b1, b2.cpu()))
                    self.assertRaises(RuntimeError, lambda: torch.bmm(b1.cpu(), b2))
                    self.assertRaises(
                        RuntimeError, lambda: torch.bmm(b1, b2, out=res2.cpu())
                    )