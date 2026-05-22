def score_tangential_damping(cond_score: torch.Tensor, uncond_score: torch.Tensor) -> torch.Tensor:
    """Drop tangential components from uncond score to align with cond score."""
    # (B, 1, ...)
    batch_num = cond_score.shape[0]
    cond_score_flat = cond_score.reshape(batch_num, 1, -1).float()
    uncond_score_flat = uncond_score.reshape(batch_num, 1, -1).float()

    # Score matrix A (B, 2, ...)
    score_matrix = torch.cat((uncond_score_flat, cond_score_flat), dim=1)
    try:
        _, _, Vh = torch.linalg.svd(score_matrix, full_matrices=False)
    except RuntimeError:
        # Fallback to CPU
        _, _, Vh = torch.linalg.svd(score_matrix.cpu(), full_matrices=False)

    # Drop the tangential components
    v1 = Vh[:, 0:1, :].to(uncond_score_flat.device)  # (B, 1, ...)
    uncond_score_td = (uncond_score_flat @ v1.transpose(-2, -1)) * v1
    return uncond_score_td.reshape_as(uncond_score).to(uncond_score.dtype)