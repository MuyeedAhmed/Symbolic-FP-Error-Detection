def _build_freqs_cis(inv_freq: torch.Tensor, attention_scaling: float, position_ids: torch.Tensor, dtype: torch.dtype,
) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
    inv_freq_e = inv_freq[None, :, None].float().expand(position_ids.shape[0], -1, 1)
    pos_e = position_ids[:, None, :].float()
    freqs = (inv_freq_e @ pos_e).transpose(1, 2)
    emb = torch.cat((freqs, freqs), dim=-1)
    cos = (emb.cos() * attention_scaling).to(dtype).unsqueeze(1)
    sin = (emb.sin() * attention_scaling).to(dtype).unsqueeze(1)
    sin_split = sin.shape[-1] // 2
    return cos, sin[..., :sin_split], -sin[..., sin_split:]