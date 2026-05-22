def compute_maxsim(q_emb: np.ndarray, d_emb: np.ndarray) -> float:
    """Compute ColBERT-style MaxSim score between query and document."""
    sim = q_emb @ d_emb.T
    return float(sim.max(axis=-1).sum())