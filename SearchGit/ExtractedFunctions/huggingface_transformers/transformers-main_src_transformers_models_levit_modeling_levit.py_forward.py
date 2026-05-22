    def forward(self, hidden_state):
        batch_size, seq_length, _ = hidden_state.shape
        key, value = (
            self.keys_values(hidden_state)
            .view(batch_size, seq_length, self.num_attention_heads, -1)
            .split([self.key_dim, self.attention_ratio * self.key_dim], dim=3)
        )
        key = key.permute(0, 2, 1, 3)
        value = value.permute(0, 2, 1, 3)

        query = self.queries(self.queries_subsample(hidden_state))
        query = query.view(batch_size, self.resolution_out**2, self.num_attention_heads, self.key_dim).permute(
            0, 2, 1, 3
        )

        attention = query @ key.transpose(-2, -1) * self.scale + self.get_attention_biases(hidden_state.device)
        attention = attention.softmax(dim=-1)
        hidden_state = (attention @ value).transpose(1, 2).reshape(batch_size, -1, self.out_dim_projection)
        hidden_state = self.projection(self.activation(hidden_state))
        return hidden_state