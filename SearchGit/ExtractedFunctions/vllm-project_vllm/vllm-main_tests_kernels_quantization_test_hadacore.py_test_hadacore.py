def test_hadacore(batch_size, hidden_dim, dtype=torch.bfloat16, device="cuda"):
    x = torch.eye(hidden_dim, dtype=dtype, device=device)
    hadamard = deterministic_hadamard_matrix(
        hidden_dim, dtype=torch.float64, device="cuda"
    ) / math.sqrt(hidden_dim)

    y = ops.hadacore_transform(x.clone())
    y_true = (x.to(hadamard.dtype) @ hadamard.T).to(y.dtype)
    assert torch.allclose(y, y_true)

    y = ops.hadacore_transform(y)
    assert torch.allclose(y, x)