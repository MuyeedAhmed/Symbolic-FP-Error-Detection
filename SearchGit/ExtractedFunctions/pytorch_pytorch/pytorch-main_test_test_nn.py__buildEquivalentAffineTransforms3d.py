def _buildEquivalentAffineTransforms3d(device, input_size, output_size, angle_rad, axis_vector):
    input_center = [(x - 1) / 2.0 for x in input_size]
    output_center = [(x - 1) / 2.0 for x in output_size]

    s = math.sin(angle_rad)
    c = math.cos(angle_rad)
    c1 = 1 - c

    intrans_ary = np.array([
        [1, 0, 0, input_center[2]],
        [0, 1, 0, input_center[3]],
        [0, 0, 1, input_center[4]],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    inscale_ary = np.array([
        [input_center[2], 0, 0, 0],
        [0, input_center[3], 0, 0],
        [0, 0, input_center[4], 0],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    l, m, n = axis_vector
    scipyRotation_ary = np.array([
        [l * l * c1 + c, m * l * c1 - n * s, n * l * c1 + m * s, 0],
        [l * m * c1 + n * s, m * m * c1 + c, n * m * c1 - l * s, 0],
        [l * n * c1 - m * s, m * n * c1 + l * s, n * n * c1 + c, 0],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    z, y, x = axis_vector
    torchRotation_ary = np.array([
        [x * x * c1 + c, y * x * c1 - z * s, z * x * c1 + y * s, 0],
        [x * y * c1 + z * s, y * y * c1 + c, z * y * c1 - x * s, 0],
        [x * z * c1 - y * s, y * z * c1 + x * s, z * z * c1 + c, 0],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    outscale_ary = np.array([
        [1.0 / output_center[2], 0, 0, 0],
        [0, 1.0 / output_center[3], 0, 0],
        [0, 0, 1.0 / output_center[4], 0],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    outtrans_ary = np.array([
        [1, 0, 0, -output_center[2]],
        [0, 1, 0, -output_center[3]],
        [0, 0, 1, -output_center[4]],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    reorder_ary = np.array([
        [0, 0, 1, 0],
        [0, 1, 0, 0],
        [1, 0, 0, 0],
        [0, 0, 0, 1],
    ], dtype=np.float64)

    transform_ary = np.dot(np.dot(np.dot(np.dot(
        intrans_ary,
        inscale_ary),
        np.linalg.inv(scipyRotation_ary)),
        outscale_ary),
        outtrans_ary)
    grid_ary = np.dot(np.dot(np.dot(reorder_ary, np.linalg.inv(scipyRotation_ary)), outscale_ary), outtrans_ary)

    transform_tensor = torch.from_numpy(torchRotation_ary).to(device, torch.float32)
    transform_tensor = transform_tensor[:3].unsqueeze(0)

    return transform_tensor, transform_ary, grid_ary