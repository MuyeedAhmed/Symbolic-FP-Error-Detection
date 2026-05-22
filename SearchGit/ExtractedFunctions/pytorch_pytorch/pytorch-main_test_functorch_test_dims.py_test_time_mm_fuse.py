    def test_time_mm_fuse(self):
        i, j, k = dims()
        A = torch.rand(3, 4)
        B = torch.rand(4, 5)

        for _ in range(10):
            r0 = A @ B

        for _ in range(10):
            a = A[i, k]
            b = B[k, j]
            r1 = (a * b).sum(k)

        with measure("pp"):
            for _ in range(10000):
                A @ B
        # magic_trace_stop_indicator()

        with measure("fc"):
            for _ in range(10000):
                (A[i, k] * B[k, j]).sum(k).order(i, j)

        with magic_trace("f.fxt"):
            for _ in range(10000):
                (A[i, k] * B[k, j]).sum(k).order(i, j)

        with magic_trace("p.fxt"):
            for _ in range(10000):
                A @ B

        # magic_trace_stop_indicator()

        torch.testing.assert_close(r1.order(i, j), r0)