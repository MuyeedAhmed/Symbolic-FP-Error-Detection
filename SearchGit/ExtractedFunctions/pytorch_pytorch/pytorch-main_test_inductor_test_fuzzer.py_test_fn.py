    def test_fn() -> bool:
        matmul_output = inp @ weight
        torch.nn.LayerNorm(hidden_size, device=GPU_TYPE)(matmul_output)
        return True