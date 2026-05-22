        def transformer(x, mask, *layer_params):
            T, C = x.shape
            params_per_layer = (
                9  # rms_w1, rms_w2, w_q, w_k, w_v, w_o, w_gate, w_up, w_down
            )

            for i in range(num_layers):
                offset = i * params_per_layer
                rms_w1 = layer_params[offset]
                rms_w2 = layer_params[offset + 1]
                w_q = layer_params[offset + 2]
                w_k = layer_params[offset + 3]
                w_v = layer_params[offset + 4]
                w_o = layer_params[offset + 5]
                w_gate = layer_params[offset + 6]
                w_up = layer_params[offset + 7]
                w_down = layer_params[offset + 8]

                # Pre-attention RMSNorm
                variance = x.pow(2).mean(-1, keepdim=True)
                h = x * torch.rsqrt(variance + 1e-6) * rms_w1

                # Multi-head self-attention
                q = (h @ w_q).view(T, num_heads, head_dim).permute(1, 0, 2)
                k = (h @ w_k).view(T, num_heads, head_dim).permute(1, 0, 2)
                v = (h @ w_v).view(T, num_heads, head_dim).permute(1, 0, 2)

                scale = 1.0 / (head_dim**0.5)
                att = (q @ k.transpose(-2, -1)) * scale
                att = att + mask
                att = torch.softmax(att, dim=-1)
                attn_out = (att @ v).permute(1, 0, 2).contiguous().view(T, C)

                x = x + (attn_out @ w_o)

                # Pre-FFN RMSNorm
                variance = x.pow(2).mean(-1, keepdim=True)
                h = x * torch.rsqrt(variance + 1e-6) * rms_w2

                # SwiGLU FFN
                gate = torch.nn.functional.silu(h @ w_gate)
                up = h @ w_up
                x = x + ((gate * up) @ w_down)

            return x