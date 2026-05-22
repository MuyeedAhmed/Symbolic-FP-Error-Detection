    def forward(
        self,
        hidden_states: torch.Tensor,
        height: int,
        width: int,
        output_attentions: bool = False,
    ) -> tuple[torch.Tensor]:
        batch_size, seq_len, num_channels = hidden_states.shape
        query_layer = self.transpose_for_scores(self.query(hidden_states))

        if self.linear_attention:
            hidden_states = hidden_states.permute(0, 2, 1).reshape(batch_size, num_channels, height, width)
            hidden_states = (
                self.spatial_reduction(self.pool(hidden_states)).reshape(batch_size, num_channels, -1).permute(0, 2, 1)
            )
            hidden_states = self.act(self.layer_norm(hidden_states))
        elif self.spatial_reduction_ratio > 1:
            hidden_states = hidden_states.permute(0, 2, 1).reshape(batch_size, num_channels, height, width)
            hidden_states = (
                self.spatial_reduction(hidden_states).reshape(batch_size, num_channels, -1).permute(0, 2, 1)
            )
            hidden_states = self.layer_norm(hidden_states)

        key_layer = self.transpose_for_scores(self.key(hidden_states))
        value_layer = self.transpose_for_scores(self.value(hidden_states))

        # Take the dot product between "query" and "key" to get the raw attention scores.
        attention_scores = torch.matmul(query_layer, key_layer.transpose(-1, -2))

        attention_scores = attention_scores / math.sqrt(self.attention_head_size)

        # Normalize the attention scores to probabilities.
        attention_probs = nn.functional.softmax(attention_scores, dim=-1)

        # This is actually dropping out entire tokens to attend to, which might
        # seem a bit unusual, but is taken from the original Transformer paper.
        attention_probs = self.attn_drop(attention_probs)
        context_layer = (attention_probs @ value_layer).transpose(1, 2).reshape(batch_size, seq_len, num_channels)
        context_layer = self.proj(context_layer)
        context_layer = self.proj_drop(context_layer)

        outputs = (context_layer, attention_probs) if output_attentions else (context_layer,)

        return outputs