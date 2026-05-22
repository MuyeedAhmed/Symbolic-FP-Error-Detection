    def forward(self, hidden_states: torch.Tensor, lm_head_weight: torch.Tensor) -> torch.Tensor:
        batch, seq_len = hidden_states.shape[:2]
        centroid_logits = self.centroids(hidden_states)

        _, top_k_indices = torch.topk(centroid_logits, k=self.centroid_intermediate_top_k, dim=-1)
        token_ordering = self.token_ordering.long()
        canonical_positions_per_cluster = token_ordering.view(self.num_centroids, self.vocab_size_per_centroid)

        # For selected top-K clusters, get canonical positions
        selected_canonical = canonical_positions_per_cluster[top_k_indices]  # [B, L, top_k, K]

        # Gather embeddings from lm_head at these canonical positions
        selected_flat = selected_canonical.reshape(-1)  # [B*L*top_k*K]
        selected_embeddings = lm_head_weight[selected_flat].view(
            batch, seq_len, self.centroid_intermediate_top_k * self.vocab_size_per_centroid, self.hidden_size
        )

        # Compute dot products: [B, L, 1, D] @ [B, L, D, top_k*K] -> [B, L, top_k*K]
        selected_logits = (hidden_states.unsqueeze(-2) @ selected_embeddings.transpose(-1, -2)).squeeze(-2)
        mask_value = selected_logits.min().item() - 1.0

        # Scatter logits directly to canonical positions in the output
        output = torch.full(
            (batch, seq_len, self.vocab_size),
            fill_value=mask_value,
            dtype=hidden_states.dtype,
            device=hidden_states.device,
        )
        scatter_idx = selected_canonical.view(batch, seq_len, -1)  # [B, L, top_k*K]
        return output.scatter_(dim=-1, index=scatter_idx, src=selected_logits)