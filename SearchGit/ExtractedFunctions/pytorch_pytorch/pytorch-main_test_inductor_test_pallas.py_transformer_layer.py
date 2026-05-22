        def transformer_layer(
            x,
            rms_w1,
            rms_w2,
            w_q,
            w_k,
            w_v,
            w_o,
            w_gate,
            w_up,
            w_down,
            mask,
        ):
            T, C = x.shape

            # Pre-attention RMSNorm
            variance = x.pow(2).mean(-1, keepdim=True)
            h = x * torch.rsqrt(variance + 1e-6) * rms_w1

            # Multi-head self-attention
            q = (h @ w_q).view(T, num_heads, head_dim).permute(1, 0, 2)  # (H, T, D)
            k = (h @ w_k).view(T, num_heads, head_dim).permute(1, 0, 2)
            v = (h @ w_v).view(T, num_heads, head_dim).permute(1, 0, 2)

            scale = 1.0 / (head_dim**0.5)
            att = (q @ k.transpose(-2, -1)) * scale  # (H, T, T)
            att = att + mask  # causal mask broadcasts (T, T) -> (H, T, T)
            att = torch.softmax(att, dim=-1)
            attn_out = (att @ v).permute(1, 0, 2).contiguous().view(T, C)  # (T, C)

            x = x + (attn_out @ w_o)

            # Pre-FFN RMSNorm
            variance = x.pow(2).mean(-1, keepdim=True)
            h = x * torch.rsqrt(variance + 1e-6) * rms_w2

            # SwiGLU FFN
            gate = torch.nn.functional.silu(h @ w_gate)
            up = h @ w_up
            x = x + ((gate * up) @ w_down)

            return x