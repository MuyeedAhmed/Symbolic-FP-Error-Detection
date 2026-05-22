        def sfdp_pattern_9_v2(query, key, value, training):
            q = query.permute(0, 2, 1, 3)
            k = key.permute(0, 2, 1, 3)
            v = value.permute(0, 2, 1, 3)
            q = q / (math.sqrt(q.size(-1)) + 0.1)
            div = q @ k.transpose(-2, -1)
            div = div.to(torch.float32)
            attn_weight = torch.softmax(div, dim=-1)
            # very low dropout to make test pass
            attn_weight = torch.dropout(attn_weight, 0.00000000001, training)
            attn_weight = attn_weight.to(torch.float16)
            return attn_weight @ v