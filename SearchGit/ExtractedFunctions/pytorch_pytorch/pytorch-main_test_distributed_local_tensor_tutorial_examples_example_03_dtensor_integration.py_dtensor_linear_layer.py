def dtensor_linear_layer(world_size: int = 4):
    """Simulate a distributed linear layer forward pass.

    Returns: (actual, expected)
    """
    batch_size, in_features, out_features = 16, 8, 4

    with LocalTensorMode(world_size):
        mesh = init_device_mesh("cpu", (world_size,))

        x = torch.randn(batch_size, in_features)
        w = torch.randn(in_features, out_features)
        b = torch.randn(out_features)

        dx = distribute_tensor(x, mesh, [Shard(0)])
        dw = distribute_tensor(w, mesh, [Replicate()])
        db = distribute_tensor(b, mesh, [Replicate()])

        dy = torch.relu(dx @ dw + db)

        expected = torch.relu(x @ w + b)
        actual = dy.full_tensor().reconcile()

    return actual, expected