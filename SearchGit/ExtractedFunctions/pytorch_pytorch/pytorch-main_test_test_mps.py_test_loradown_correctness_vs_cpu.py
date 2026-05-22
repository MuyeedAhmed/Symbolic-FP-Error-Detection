    def test_loradown_correctness_vs_cpu(self, padding, vector_dim):
        torch.manual_seed(13)

        base_size = 64
        physical_size = base_size + padding

        a_mps = torch.rand(vector_dim, base_size, device='mps', dtype=torch.half)
        b_mps = torch.rand(vector_dim, physical_size, device='mps', dtype=torch.half)[:, :base_size].t()

        a_cpu = a_mps.cpu()
        b_cpu = b_mps.cpu()

        result_cpu = (a_cpu @ b_cpu)
        result_mps = (a_mps @ b_mps).cpu()

        torch.testing.assert_close(result_mps, result_cpu, rtol=1e-3, atol=1e-3)