    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        batch_size, sequence_length, _ = inputs.shape
        projection = self.dense(inputs).reshape(batch_size, sequence_length, 2, self.head_size)
        projection = self.dropout(projection)
        queries, keys = torch.unbind(projection, dim=2)
        logits = (queries @ keys.transpose(-2, -1)) / (self.head_size**0.5)
        mask = torch.tril(torch.ones(sequence_length, sequence_length, device=logits.device)).bool()
        return logits.masked_fill(mask.unsqueeze(0), -1e4)