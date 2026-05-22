    def forward(self, hidden_state, output_attentions=False):
        batch_size, height, width, _ = hidden_state.shape
        # qkv with shape (3, batch_size, num_heads, height * width, num_channels)
        qkv = self.qkv(hidden_state).reshape(batch_size, height * width, 3, self.num_heads, -1).permute(2, 0, 3, 1, 4)
        # queries, keys and values have shape (batch_size * num_heads, height * width, num_channels)
        queries, keys, values = qkv.reshape(3, batch_size * self.num_heads, height * width, -1).unbind(0)

        attention_scores = (queries * self.scale) @ keys.transpose(-2, -1)

        if self.use_relative_position_embeddings:
            attention_scores = add_decomposed_relative_positions(
                attention_scores, queries, self.rel_pos_h, self.rel_pos_w, (height, width), (height, width)
            )

        attention_probs = attention_scores.softmax(dim=-1)

        hidden_state = attention_probs @ values
        hidden_state = hidden_state.view(batch_size, self.num_heads, height, width, -1)
        hidden_state = hidden_state.permute(0, 2, 3, 1, 4)
        hidden_state = hidden_state.reshape(batch_size, height, width, -1)
        hidden_state = self.proj(hidden_state)

        if output_attentions:
            attention_probs = attention_probs.reshape(
                batch_size, self.num_heads, attention_probs.shape[-2], attention_probs.shape[-1]
            )
            outputs = (hidden_state, attention_probs)
        else:
            outputs = (hidden_state,)

        return outputs