def precompute_partial_rope(head_dim, rotary_dim, position_ids, theta, device=None, mrope_section=None):
    """Compute RoPE frequencies for partial rotary embeddings."""
    theta_numerator = torch.arange(0, rotary_dim, 2, device=device).float()
    inv_freq = 1.0 / (theta ** (theta_numerator / rotary_dim))

    inv_freq_expanded = inv_freq[None, :, None].float().expand(position_ids.shape[0], -1, 1)
    position_ids_expanded = position_ids[:, None, :].float()
    freqs = (inv_freq_expanded.float() @ position_ids_expanded.float()).transpose(1, 2)
    emb = torch.cat((freqs, freqs), dim=-1)
    cos = emb.cos()
    sin = emb.sin()

    if mrope_section is not None and position_ids.shape[0] == 3:
        mrope_section_2 = [s * 2 for s in mrope_section]
        cos = torch.cat([m[i % 3] for i, m in enumerate(cos.split(mrope_section_2, dim=-1))], dim=-1).unsqueeze(0)
        sin = torch.cat([m[i % 3] for i, m in enumerate(sin.split(mrope_section_2, dim=-1))], dim=-1).unsqueeze(0)

    cos = cos.unsqueeze(1)
    sin = sin.unsqueeze(1)
    sin_split = sin.shape[-1] // 2
    return (cos, sin[..., :sin_split], -sin[..., sin_split:])