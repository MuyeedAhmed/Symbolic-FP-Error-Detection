    def forward(self, hidden_states):
        batch_size, num, channel = hidden_states.shape
        qkv = (
            self.qkv(hidden_states)
            .reshape(batch_size, num, 3, self.num_heads, channel // self.num_heads)
            .permute(2, 0, 3, 1, 4)
        )
        query, key, value = qkv[0], qkv[1], qkv[2]

        attention_probs = (query @ key.transpose(-2, -1)) * self.scale
        attention_probs = attention_probs.softmax(dim=-1)
        attention_probs = self.attn_drop(attention_probs)

        context_layer = (attention_probs @ value).transpose(1, 2).reshape(batch_size, num, channel)
        context_layer = self.proj(context_layer)
        context_layer = self.proj_drop(context_layer)
        return (context_layer, attention_probs)