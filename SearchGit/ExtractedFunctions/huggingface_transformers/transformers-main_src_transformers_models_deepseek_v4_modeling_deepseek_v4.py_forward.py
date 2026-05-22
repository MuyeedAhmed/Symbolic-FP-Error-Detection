    def forward(self, x, position_ids, layer_type=None):
        # Key difference vs Laguna's forward: no `torch.cat([freqs, freqs], dim=-1)`
        # duplication. V4's interleaved RoPE pairs consecutive channels, so we only need
        # `rope_head_dim // 2` unique θ entries — the `apply_rotary_pos_emb` helper does
        # the `repeat_interleave(2)` next to the rotation math, where the link between
        # the doubled dim and `rotate_half` is local and obvious.
        inv_freq = getattr(self, f"{layer_type}_inv_freq")
        attention_scaling = getattr(self, f"{layer_type}_attention_scaling")
        inv_freq_expanded = inv_freq[None, :, None].float().expand(position_ids.shape[0], -1, 1).to(x.device)
        position_ids_expanded = position_ids[:, None, :].float()
        device_type = x.device.type if isinstance(x.device.type, str) and x.device.type != "mps" else "cpu"
        with maybe_autocast(device_type=device_type, enabled=False):
            freqs = (inv_freq_expanded.float() @ position_ids_expanded.float()).transpose(1, 2)
            cos = freqs.cos() * attention_scaling
            sin = freqs.sin() * attention_scaling
        return cos.to(dtype=x.dtype), sin.to(dtype=x.dtype)