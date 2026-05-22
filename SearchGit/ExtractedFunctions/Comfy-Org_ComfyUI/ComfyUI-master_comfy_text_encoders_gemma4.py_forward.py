    def forward(self, x, position_embeddings=None, attn_mask=None):
        B, S, _ = x.shape

        q = self.q_proj(x).float().view(B, S, self.num_heads, self.head_dim)
        k = self.k_proj(x).float().view(B, S, self.num_heads, self.head_dim)
        v = self.v_proj(x).float().view(B, S, self.num_heads, self.head_dim)

        q = q * self.q_scale * torch.nn.functional.softplus(self.per_dim_scale)
        k = k * self.k_scale

        q_blocks = self._convert_to_block(q)
        k_context = self._extract_block_context(k)
        v_context = self._extract_block_context(v)
        num_blocks = q_blocks.shape[1]

        rel_k = self.relative_k_proj(position_embeddings).view(-1, self.num_heads, self.head_dim).to(q.dtype)

        queries = q_blocks.permute(0, 3, 1, 2, 4)  # [B, H, NB, CS, D]
        matrix_ac = queries @ k_context.permute(0, 3, 1, 4, 2)

        queries_flat = queries.reshape(B, self.num_heads, -1, self.head_dim)
        matrix_bd = queries_flat @ rel_k.permute(1, 2, 0)
        matrix_bd = matrix_bd.reshape(B, self.num_heads, num_blocks, self.chunk_size, -1)
        matrix_bd = self._rel_shift(matrix_bd)

        attn_weights = matrix_ac + matrix_bd
        attn_weights = torch.tanh(attn_weights / self.softcap) * self.softcap

        # Mask out invalid positions in chunk context (matching reference's masked_fill approach)
        if attn_mask is None:
            attn_mask = self._build_blocked_mask(S, num_blocks, x.device)
        attn_weights = attn_weights.masked_fill(attn_mask.logical_not(), -1e9)

        attn_weights = torch.nn.functional.softmax(attn_weights, dim=-1, dtype=torch.float32).to(v.dtype)
        out = attn_weights @ v_context.permute(0, 3, 1, 2, 4)
        out = out.permute(0, 2, 3, 1, 4).reshape(B, num_blocks * self.chunk_size, -1)
        out = out[:, :S].contiguous()
        return self.post(out.to(self.post.linear.weight.dtype))