        def forward(
            x,
            q_proj,
            k_proj,
            v_proj,
            o_proj,
            embed_norm,
            hidden_norm,
            cos,
            sin,
        ):
            batch, seq_len, _ = x.shape

            # Eagle3 first layer: split concatenated [embeds, hidden] input
            mid = x.shape[2] // 2
            embeds, hidden = x.split(mid, dim=-1)

            # Dual RMSNorm (pow, sum, div, mul in backward)
            embeds = embed_norm(embeds)
            hidden = hidden_norm(hidden)
            residual = hidden

            # Recombine for attention input (2 * HIDDEN_SIZE)
            x = torch.cat([embeds, hidden], dim=-1)

            # Adding a graph break here "fixes" the issue
            # by breaking up the fused op
            # torch._dynamo.graph_break()

            # Q/K/V projections from 2*hidden_size input
            q = q_proj(x).view(batch, seq_len, NUM_HEADS, HEAD_DIM).transpose(1, 2)
            k = k_proj(x).view(batch, seq_len, NUM_KV_HEADS, HEAD_DIM).transpose(1, 2)
            v = v_proj(x).view(batch, seq_len, NUM_KV_HEADS, HEAD_DIM).transpose(1, 2)

            q, k = apply_rotary_pos_emb(q, k, cos, sin)
            k = torch.repeat_interleave(k, NUM_HEADS // NUM_KV_HEADS, dim=1)
            v = torch.repeat_interleave(v, NUM_HEADS // NUM_KV_HEADS, dim=1)
            out = q.contiguous() @ k.contiguous().transpose(-2, -1) @ v.contiguous()

            out = out.transpose(1, 2).contiguous().reshape(batch, seq_len, -1)
            return o_proj(out) + residual