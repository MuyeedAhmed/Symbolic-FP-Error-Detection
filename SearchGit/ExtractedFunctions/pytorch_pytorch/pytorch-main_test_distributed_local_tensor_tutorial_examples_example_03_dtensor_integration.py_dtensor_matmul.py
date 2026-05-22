def dtensor_matmul(world_size: int = 4):
    """Perform matrix multiplication with DTensors.

    Returns: (actual, expected)
    """
    with LocalTensorMode(world_size):
        mesh = init_device_mesh("cpu", (world_size,))

        a = torch.randn(8, 4)
        b = torch.randn(4, 6)

        da = distribute_tensor(a, mesh, [Shard(0)])
        db = distribute_tensor(b, mesh, [Replicate()])

        dc = da @ db

        expected = a @ b
        actual = dc.full_tensor().reconcile()

    return actual, expected