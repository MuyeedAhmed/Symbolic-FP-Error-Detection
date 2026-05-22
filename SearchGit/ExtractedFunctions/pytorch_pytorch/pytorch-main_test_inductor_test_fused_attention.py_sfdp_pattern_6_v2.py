        def sfdp_pattern_6_v2(query, key, value, training):
            attn_mask = torch.ones(
                query.size(-2), key.size(-2), dtype=torch.bool, device=query.device
            ).tril(diagonal=0)
            attn_mask = attn_mask.masked_fill(
                torch.logical_not(attn_mask), -float("inf")
            )
            attn_weight = torch.softmax(
                (query @ key.transpose(-2, -1) / (math.sqrt(query.size(-1)) + 0.1))
                + attn_mask,
                dim=-1,
            )
            attn_weight = torch.nn.functional.dropout(attn_weight, 0.5, training)
            return attn_weight @ value