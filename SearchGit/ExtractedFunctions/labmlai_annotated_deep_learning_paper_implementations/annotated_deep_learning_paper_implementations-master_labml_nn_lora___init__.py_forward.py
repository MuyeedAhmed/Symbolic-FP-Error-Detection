    def forward(self, x: torch.Tensor):
        # Compute the embeddings $\text{onehot}(x) W_0$
        result = nn.functional.embedding(x, self.weight)

        # Add $\frac{\alpha}{r} \text{onehot}(x) \Delta W^T = \frac{\alpha}{r} \text{onehot}(x) A^T B^T$
        result += (nn.functional.embedding(x, self.lora_a.T) @ self.lora_b.T) * self.scaling

        #
        return result