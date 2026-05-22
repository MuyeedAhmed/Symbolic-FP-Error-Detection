    def test_tau_token_output_range(self):
        """Test tau token scaling output is bounded by tanh."""
        import torch.nn.functional as F

        seq_len = 100
        qkv_dim = 6144  # 2048 * 3
        num_heads = 32

        qkv = torch.randn(seq_len, qkv_dim)
        tau_wq = torch.randn(num_heads, qkv_dim)

        tok_feat = F.gelu(qkv)
        tok_q = torch.tanh(tok_feat @ tau_wq.t())

        assert tok_q.shape == (seq_len, num_heads)
        # tanh output is bounded by [-1, 1]
        assert tok_q.min() >= -1.0
        assert tok_q.max() <= 1.0