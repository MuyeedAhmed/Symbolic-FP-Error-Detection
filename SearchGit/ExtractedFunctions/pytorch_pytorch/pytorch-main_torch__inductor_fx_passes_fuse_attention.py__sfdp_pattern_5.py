def _sfdp_pattern_5(query, key, value, attn_mask, inv_scale):
    attn_weight = torch.softmax(
        (query @ key.transpose(-2, -1) / (inv_scale)) + attn_mask, dim=-1
    )
    # attn_weight = torch.dropout(attn_weight, dropout_p)
    return attn_weight @ value