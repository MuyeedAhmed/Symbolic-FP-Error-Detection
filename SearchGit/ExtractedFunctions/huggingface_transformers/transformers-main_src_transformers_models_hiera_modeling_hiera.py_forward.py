    def forward(
        self,
        hidden_states: torch.Tensor,
        output_attentions: bool = False,
    ) -> tuple[torch.Tensor, torch.Tensor | None]:
        """Input should be of shape [batch, tokens, channels]."""
        batch_size, seq_len, _ = hidden_states.shape

        num_windows = 1
        if self.use_mask_unit_attn:
            num_windows = seq_len // (self.query_stride * self.window_size)

        qkv = self.qkv(hidden_states)
        qkv = qkv.reshape(batch_size, -1, num_windows, 3, self.num_heads, self.head_dim)
        qkv = qkv.permute(3, 0, 4, 2, 1, 5)

        query, key, value = qkv.unbind(0)

        if self.query_stride > 1:
            # Refer to unroll to see how this performs a maxpool-Nd
            query = query.view(batch_size, self.num_heads, num_windows, self.query_stride, -1, self.head_dim)
            query = query.max(dim=3).values

        attn_weights = (query * self.scale) @ key.transpose(-1, -2)
        attn_weights = attn_weights.softmax(dim=-1)

        attn_output = attn_weights @ value
        attn_output = attn_output.transpose(1, 3).reshape(batch_size, -1, self.hidden_size_output)
        attn_output = self.proj(attn_output)

        return (attn_output, attn_weights) if output_attentions else (attn_output, None)