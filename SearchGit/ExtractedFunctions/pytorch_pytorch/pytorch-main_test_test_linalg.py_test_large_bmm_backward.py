    def test_large_bmm_backward(self, device):
        A = torch.randn([1024, 2, 1024], device="cuda").mT.contiguous().mT
        B = torch.randn([1, 1024, 65536], device="cuda", requires_grad=True)
        G = torch.randn([1024, 2, 65536], device="cuda")

        # Should not create an intermediary tensor of size [1024, 1024, 65536] (256GB of memory) and OOM
        (A @ B).backward(G)