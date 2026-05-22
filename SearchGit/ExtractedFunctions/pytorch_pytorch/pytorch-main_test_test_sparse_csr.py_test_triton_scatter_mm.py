    def test_triton_scatter_mm(self, device, dtype):
        from torch.sparse._triton_ops import scatter_mm
        from functools import partial
        tensor = partial(make_tensor, device=device, dtype=dtype, low=0.5, high=1.5)
        sizes = [8, 16]
        for m, k, n in itertools.product(sizes, sizes, sizes):
            blocks = torch.stack([tensor(m, k), tensor(m, k)])
            others = torch.stack([tensor(k, n), tensor(k, n)])

            expected = torch.stack([blocks[0] @ others[0] + blocks[1] @ others[0],
                                    blocks[0] @ others[1],
                                    blocks[1] @ others[1]])

            indices_data = (
                'scatter_mm',
                torch.tensor([0, 2, 3, 4], dtype=torch.int32, device=device),
                torch.tensor([[0, 0], [1, 0], [0, 1], [1, 1]], dtype=torch.int32, device=device))

            result = scatter_mm(blocks, others, indices_data=indices_data)

            self.assertEqual(result, expected)

            indices_data = (
                'bsr_strided_mm',
                torch.tensor([0, 2, 4, 5, 6], dtype=torch.int32, device=device),
                torch.tensor([0, n, 2 * n * m, 2 * n * m + n], dtype=torch.int32, device=device),
                torch.tensor([1, 0, 1, 0, 1, 1], dtype=torch.int32, device=device),
                torch.tensor([0, 2 * k * n, n, 2 * k * n + n, 2 * k * n, 2 * k * n + n],
                             dtype=torch.int32, device=device),
                dict(SPLIT_N=2, is_compressed=False, TILE_M=m, TILE_N=n, GROUP_SIZE=1)
            )

            for bsize in [(), (2,), (3, 4)]:
                other = tensor(*bsize, 2 * k, 2 * n)
                expected = torch.cat([
                    torch.cat([blocks[1], blocks[0]], dim=1),
                    torch.cat([torch.zeros_like(blocks[0]), blocks[1]], dim=1)], dim=0) @ other
                result = scatter_mm(blocks, other, indices_data=indices_data)
                self.assertEqual(result, expected)