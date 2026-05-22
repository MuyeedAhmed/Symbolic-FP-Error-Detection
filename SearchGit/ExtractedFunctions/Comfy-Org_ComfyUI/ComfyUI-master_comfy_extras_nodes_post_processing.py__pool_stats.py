    def _pool_stats(images, device, is_reinhard, eps):
        """Two-pass pooled mean + std/cov across all frames."""
        N, C = images.shape[0], images.shape[3]
        HW = images.shape[1] * images.shape[2]
        mean = torch.zeros(C, 1, device=device, dtype=torch.float32)
        for i in range(N):
            mean += ColorTransfer._to_lab(images, i, device).view(C, -1).mean(dim=-1, keepdim=True)
        mean /= N
        acc = torch.zeros(C, 1 if is_reinhard else C, device=device, dtype=torch.float32)
        for i in range(N):
            centered = ColorTransfer._to_lab(images, i, device).view(C, -1) - mean
            if is_reinhard:
                acc += (centered * centered).mean(dim=-1, keepdim=True)
            else:
                acc += centered @ centered.T / HW
        if is_reinhard:
            return mean, torch.sqrt(acc / N).clamp_min_(eps)
        return mean, acc / N