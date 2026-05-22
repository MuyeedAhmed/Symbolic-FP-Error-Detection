        def sfdp_pattern_7(query, key, value, training):
            q = query.permute(0, 2, 1, 3)
            k = key.permute(0, 2, 1, 3)
            v = value.permute(0, 2, 1, 3)
            div = q @ k.transpose(-2, -1) / math.sqrt(q.size(-1))
            div = div.to(torch.float32)
            attn_weight = torch.softmax(div, dim=-1)
            # Set to False
            attn_weight = torch.dropout(attn_weight, 0.00000000001, training)
            attn_weight = attn_weight.to(torch.float16)
            return attn_weight @ v