def hybrid_parallelism():
    """Combine data parallel and tensor parallel.

    Returns: (actual, expected)
    """
    world_size = 8
    dp_size, tp_size = 4, 2

    with LocalTensorMode(world_size):
        mesh = init_device_mesh("cpu", (dp_size, tp_size), mesh_dim_names=("dp", "tp"))

        x = torch.randn(16, 8)
        dx = distribute_tensor(x, mesh, [Shard(0), Replicate()])

        w = torch.randn(8, 12)
        dw = distribute_tensor(w, mesh, [Replicate(), Shard(1)])

        dy = dx @ dw

        expected = x @ w
        actual = dy.full_tensor().reconcile()

    return actual, expected