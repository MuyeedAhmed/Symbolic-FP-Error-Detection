    def forward(self, hidden_states: torch.Tensor, router_indices=None, routing_weights=None) -> torch.Tensor:
        original_dtype = hidden_states.dtype

        # Accumulate over fp32
        next_states = torch.zeros_like(hidden_states, dtype=torch.float32, device=hidden_states.device)
        with torch.no_grad():
            expert_mask = torch.nn.functional.one_hot(
                router_indices, num_classes=self.num_experts
            )  # masking is also a class
            expert_mask = expert_mask.permute(2, 1, 0)
            # we sum on the top_k and on the sequence length to get which experts
            # are hit this time around
            expert_hit = torch.greater(expert_mask.sum(dim=(-1, -2)), 0).nonzero()

        # Key change to original gpt oss is to stay in fp32 precision for all linear projections / muls
        for expert_idx in expert_hit:
            # expert_idx only have 1 element, so we can use scale for fast indexing
            expert_idx = expert_idx[0]
            # skip masking index
            if expert_idx == self.num_experts:
                continue
            top_k_pos, token_idx = torch.where(expert_mask[expert_idx])
            current_state = hidden_states[token_idx]
            gate_up = (
                current_state.float() @ self.gate_up_proj[expert_idx].float()
                + self.gate_up_proj_bias[expert_idx].float()
            )
            gated_output = self._apply_gate(gate_up).float()
            out = gated_output.float() @ self.down_proj[expert_idx].float() + self.down_proj_bias[expert_idx].float()
            weighted_output = out * routing_weights[token_idx, top_k_pos, None].float()
            next_states.index_add_(0, token_idx, weighted_output)

        return next_states.to(original_dtype)