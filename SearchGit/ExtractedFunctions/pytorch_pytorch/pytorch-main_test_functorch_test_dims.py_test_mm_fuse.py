    def test_mm_fuse(self):
        i, j, k = dims()
        A = torch.rand(3, 4)
        B = torch.rand(4, 5)

        C = (A[i, k] * B[k, j]).sum(k).order(i, j)
        torch.testing.assert_close(C, A @ B)