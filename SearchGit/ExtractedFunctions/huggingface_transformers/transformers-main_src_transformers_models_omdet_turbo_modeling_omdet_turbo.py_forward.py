    def forward(self, hidden_states, mask=None, encode_type="task"):
        text_outputs = self.model(hidden_states)
        pooled_output = text_outputs[0]
        if encode_type == "task":
            if mask is None:
                raise ValueError("mask is required for task encoding")
            max_len = (mask != 0).sum(1).max().item()
            truncated_mask = mask[:, :max_len]
            truncated_output = pooled_output[:, :max_len, :]
            return truncated_output.transpose(0, 1), truncated_mask
        elif encode_type == "class":
            max_pooled_output = pooled_output[torch.arange(pooled_output.shape[0]), hidden_states.argmax(dim=-1)]
            projected_output = max_pooled_output @ self.text_projection
            return projected_output
        else:
            raise ValueError(f"encode_type {encode_type} is not supported")