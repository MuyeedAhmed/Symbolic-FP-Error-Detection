    def forward(self, hidden_states, output_attentions: bool = False):
        batch_size, hidden_size, num_channels = hidden_states.shape
        qkv = (
            self.qkv(hidden_states)
            .reshape(batch_size, hidden_size, 3, self.num_heads, num_channels // self.num_heads)
            .permute(2, 0, 3, 1, 4)
        )
        query, key, value = qkv[0], qkv[1], qkv[2]

        attention_probs = (query @ key.transpose(-2, -1)) * self.scale
        attention_probs = attention_probs.softmax(dim=-1)
        attention_probs = self.attn_drop(attention_probs)

        context_layer = (attention_probs @ value).transpose(1, 2).reshape(batch_size, hidden_size, num_channels)

        outputs = (context_layer, attention_probs) if output_attentions else (context_layer,)

        return outputs