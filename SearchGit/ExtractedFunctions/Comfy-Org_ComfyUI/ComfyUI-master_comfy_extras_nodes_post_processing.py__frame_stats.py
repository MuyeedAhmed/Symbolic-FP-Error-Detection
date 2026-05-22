    def _frame_stats(lab_flat, hw, is_reinhard, eps):
        """Per-frame mean + std/cov."""
        mean = lab_flat.mean(dim=-1, keepdim=True)
        if is_reinhard:
            return mean, lab_flat.std(dim=-1, keepdim=True, unbiased=False).clamp_min_(eps)
        centered = lab_flat - mean
        return mean, centered @ centered.T / hw