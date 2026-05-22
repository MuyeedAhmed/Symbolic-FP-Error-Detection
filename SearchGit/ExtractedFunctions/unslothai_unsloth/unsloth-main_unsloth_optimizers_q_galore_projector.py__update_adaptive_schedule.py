    def _update_adaptive_schedule(
        self,
        float_ortho: torch.Tensor,
        side: str,
    ) -> None:
        """Track subspace stability and increase ``update_proj_gap`` if stable."""
        self.svd_count += 1

        if side == "right":
            current_vector = float_ortho[:1, :].flatten()
        else:
            current_vector = float_ortho[:, :1].flatten()

        if self.past_ortho_vector is not None:
            cos_sim = torch.dot(self.past_ortho_vector, current_vector).item()

            self.queue.append(cos_sim)

            if (
                len(self.queue) == self.queue.maxlen
                and sum(self.queue) / len(self.queue) >= self.cos_threshold
            ):
                self.update_proj_gap = int(self.update_proj_gap * self.gamma_proj)

        self.past_ortho_vector = current_vector.clone()