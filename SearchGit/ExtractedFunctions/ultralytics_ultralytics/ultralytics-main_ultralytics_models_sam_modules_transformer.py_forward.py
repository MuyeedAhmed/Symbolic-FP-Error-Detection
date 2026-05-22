    def forward(self, q: torch.Tensor, k: torch.Tensor, v: torch.Tensor) -> torch.Tensor:
        """Apply multi-head attention to query, key, and value tensors with optional downsampling.

        Args:
            q (torch.Tensor): Query tensor with shape (B, N_q, embedding_dim).
            k (torch.Tensor): Key tensor with shape (B, N_k, kv_in_dim).
            v (torch.Tensor): Value tensor with shape (B, N_k, kv_in_dim).

        Returns:
            (torch.Tensor): Output tensor after attention with shape (B, N_q, embedding_dim).
        """
        # Input projections
        q = self.q_proj(q)
        k = self.k_proj(k)
        v = self.v_proj(v)

        # Separate into heads
        q = self._separate_heads(q, self.num_heads)
        k = self._separate_heads(k, self.num_heads)
        v = self._separate_heads(v, self.num_heads)

        # Attention
        _, _, _, c_per_head = q.shape
        attn = q @ k.permute(0, 1, 3, 2)  # B x N_heads x N_tokens x N_tokens
        attn = attn / math.sqrt(c_per_head)
        attn = torch.softmax(attn, dim=-1)

        # Get output
        out = attn @ v
        out = self._recombine_heads(out)
        return self.out_proj(out)