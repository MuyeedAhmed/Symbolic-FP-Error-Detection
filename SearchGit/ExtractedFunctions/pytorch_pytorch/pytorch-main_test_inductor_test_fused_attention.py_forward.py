            def forward(self, query, key, value, attn_mask) -> torch.Tensor:
                attn_weight = torch.softmax(
                    query @ key.transpose(-2, -1) / math.sqrt(query.size(-1))
                    + attn_mask,
                    dim=-1,
                )
                return attn_weight @ value