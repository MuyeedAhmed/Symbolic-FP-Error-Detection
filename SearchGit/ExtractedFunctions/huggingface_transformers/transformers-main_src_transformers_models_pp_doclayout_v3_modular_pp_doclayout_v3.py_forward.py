    def forward(self, inputs):
        batch_size, sequence_length, _ = inputs.shape
        query_key_projection = self.dense(inputs).reshape(batch_size, sequence_length, 2, self.head_size)
        query_key_projection = self.dropout(query_key_projection)
        queries, keys = torch.unbind(query_key_projection, dim=2)

        logits = (queries @ keys.transpose(-2, -1)) / (self.head_size**0.5)
        mask = torch.tril(torch.ones(sequence_length, sequence_length, device=logits.device)).bool()
        logits = logits.masked_fill(mask.unsqueeze(0), -1e4)

        return logits