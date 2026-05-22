    def decode_latents(self, hidden_states):
        batch_size, hidden_dim, sequence_length = hidden_states.shape
        encodings = hidden_states.permute(0, 2, 1).reshape(batch_size * sequence_length, hidden_dim)
        codebook = self.codebook.weight  # codebook: (N x D)

        # L2 normalize encodings and codebook (ViT-VQGAN)
        encodings = F.normalize(encodings)
        codebook = F.normalize(codebook)

        # Compute euclidean distance with codebook
        l2_norm = encodings.pow(2).sum(1, keepdim=True)
        dist = -(l2_norm - 2 * encodings @ codebook.t()) + codebook.pow(2).sum(1, keepdim=True).t()

        indices = dist.max(1)[1]
        indices = indices.reshape(hidden_states.size(0), -1)
        quantized_representation = self.codebook(indices).transpose(1, 2)
        return quantized_representation, indices