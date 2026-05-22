        def transformer_block(x, w_q, w_k, w_v, w_proj, w_fc, w_out, mask):
            T, C = x.shape

            # === Self-Attention ===
            q = x @ w_q  # (T, C)
            k = x @ w_k
            v = x @ w_v

            # Scaled dot-product attention
            scale = 1.0 / (C**0.5)
            att = (q @ k.t()) * scale  # (T, T)
            att = att + mask  # Apply causal mask
            att = torch.softmax(att, dim=-1)
            attn_out = att @ v  # (T, C)

            # Output projection + residual
            x = x + (attn_out @ w_proj)

            # === MLP (Feed-Forward) ===
            h = x @ w_fc  # (T, 4C)
            # GELU activation (tanh approximation)
            h = 0.5 * h * (1.0 + torch.tanh(0.7978845608 * (h + 0.044715 * h * h * h)))
            x = x + (h @ w_out)  # Residual + project back

            return x