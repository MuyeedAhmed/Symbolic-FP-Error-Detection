    def forward(
        self,
        positions: torch.Tensor,
        hidden_states: torch.Tensor,
    ) -> torch.Tensor:
        qkv, _ = self.qkv_proj(hidden_states)

        q, k, v = qkv.split(
            [
                self.num_heads_per_partition * self.head_dim,
                self.num_kv_heads_per_partition * self.head_dim,
                self.num_kv_heads_per_partition * self.head_dim,
            ],
            dim=-1,
        )

        # Apply tau scaling to Q and V
        # Tau scaling has two components:
        # 1. Token-based: tok_q = tanh(gelu(qkv) @ tau_wq.T)
        # 2. Position-based: tau_pos = 1 + (sigmoid(alpha * log(pos+1)) - 0.5)
        # Final: tau = tok + tau_pos
        #
        # For TP, tau weights are sharded by head, but qkv_dim is kept full

        # Get full qkv for tau computation
        # With TP, reconstruct qkv in correct layout [q_full, k_full, v_full]
        # (all-gather would produce [q_0, k_0, v_0, q_1, k_1, v_1] - wrong)
        if self.tp_size > 1:
            # All-gather once, then reconstruct [q_full, k_full, v_full].
            qkv_full_sharded = tensor_model_parallel_all_gather(qkv.contiguous())
            q_local_dim = q.shape[-1]
            kv_local_dim = k.shape[-1]
            qkv_full_sharded = qkv_full_sharded.view(
                qkv.shape[0],
                self.tp_size,
                q_local_dim + 2 * kv_local_dim,
            )
            q_full = qkv_full_sharded[:, :, :q_local_dim].reshape(qkv.shape[0], -1)
            k_full = qkv_full_sharded[
                :, :, q_local_dim : q_local_dim + kv_local_dim
            ].reshape(qkv.shape[0], -1)
            v_full = qkv_full_sharded[:, :, q_local_dim + kv_local_dim :].reshape(
                qkv.shape[0], -1
            )
            qkv_full = torch.cat([q_full, k_full, v_full], dim=-1).contiguous()
        else:
            qkv_full = qkv

        # Compute tau scaling factors matching HF implementation exactly:
        # tok_feat = gelu(qkv)
        # tok_q = tanh(tok_feat @ tau_wq.T)  # [num_tokens, num_heads]
        # tau_pos = 1 + (sigmoid(alpha * log(pos+1)) - 0.5)  # [num_heads, num_tokens]
        # tau = (tok_q.T + tau_pos).T  # [num_tokens, num_heads]
        num_tokens = qkv_full.shape[0]
        orig_dtype = q.dtype

        # Token-based component
        tok_feat = F.gelu(qkv_full)  # Apply GELU activation
        tok_q = torch.tanh(tok_feat @ self.tau_wq.t())  # [N, H_per_partition]
        tok_v = torch.tanh(tok_feat @ self.tau_wv.t())  # [N, H_per_partition]

        # Position-based component
        # tau_pos = 1 + (sigmoid(alpha * log(pos+1)) - 0.5)
        # positions is [num_tokens], need to compute for each head
        # tau_alpha: [num_heads_per_partition]
        pos_float = (positions.to(orig_dtype) + 1.0).clamp(min=1e-6)
        pos_log = pos_float.log()  # [num_tokens]
        # alpha[:, None] * pos_log[None, :] -> [num_heads, num_tokens]
        tau_pos = 1.0 + (
            torch.sigmoid(self.tau_alpha[:, None] * pos_log[None, :]) - 0.5
        )  # [H_per_partition, N]

        # Combine token and position components
        tau_q = (tok_q + tau_pos.t()).to(orig_dtype)  # [N, H_per_partition]
        tau_v = (tok_v + tau_pos.t()).to(orig_dtype)  # [N, H_per_partition]

        # Reshape q and v to apply per-head tau scaling
        q = q.view(num_tokens, self.num_heads_per_partition, self.head_dim)
        v = v.view(num_tokens, self.num_kv_heads_per_partition, self.head_dim)

        # Apply tau scaling
        q = q * tau_q.unsqueeze(-1)
        v = v * tau_v[:, : self.num_kv_heads_per_partition].unsqueeze(-1)

        # Reshape back
        q = q.view(num_tokens, -1)
        v = v.view(num_tokens, -1)

        q, k = self.rotary_emb(positions, q, k)

        attn_output = self.attn(q, k, v)

        output, _ = self.out_proj(attn_output)
        return output