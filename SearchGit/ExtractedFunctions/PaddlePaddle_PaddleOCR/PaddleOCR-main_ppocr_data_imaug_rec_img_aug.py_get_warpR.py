def get_warpR(config):
    """
    get_warpR
    """
    anglex, angley, anglez, fov, w, h, r = (
        config.anglex,
        config.angley,
        config.anglez,
        config.fov,
        config.w,
        config.h,
        config.r,
    )
    if w > 69 and w < 112:
        anglex = anglex * 1.5

    z = np.sqrt(w**2 + h**2) / 2 / np.tan(rad(fov / 2))
    # Homogeneous coordinate transformation matrix
    rx = np.array(
        [
            [1, 0, 0, 0],
            [0, np.cos(rad(anglex)), -np.sin(rad(anglex)), 0],
            [
                0,
                -np.sin(rad(anglex)),
                np.cos(rad(anglex)),
                0,
            ],
            [0, 0, 0, 1],
        ],
        np.float32,
    )
    ry = np.array(
        [
            [np.cos(rad(angley)), 0, np.sin(rad(angley)), 0],
            [0, 1, 0, 0],
            [
                -np.sin(rad(angley)),
                0,
                np.cos(rad(angley)),
                0,
            ],
            [0, 0, 0, 1],
        ],
        np.float32,
    )
    rz = np.array(
        [
            [np.cos(rad(anglez)), np.sin(rad(anglez)), 0, 0],
            [-np.sin(rad(anglez)), np.cos(rad(anglez)), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
        ],
        np.float32,
    )
    r = rx.dot(ry).dot(rz)
    # generate 4 points
    pcenter = np.array([h / 2, w / 2, 0, 0], np.float32)
    p1 = np.array([0, 0, 0, 0], np.float32) - pcenter
    p2 = np.array([w, 0, 0, 0], np.float32) - pcenter
    p3 = np.array([0, h, 0, 0], np.float32) - pcenter
    p4 = np.array([w, h, 0, 0], np.float32) - pcenter
    dst1 = r.dot(p1)
    dst2 = r.dot(p2)
    dst3 = r.dot(p3)
    dst4 = r.dot(p4)
    list_dst = np.array([dst1, dst2, dst3, dst4])
    org = np.array([[0, 0], [w, 0], [0, h], [w, h]], np.float32)
    dst = np.zeros((4, 2), np.float32)
    # Project onto the image plane
    dst[:, 0] = list_dst[:, 0] * z / (z - list_dst[:, 2]) + pcenter[0]
    dst[:, 1] = list_dst[:, 1] * z / (z - list_dst[:, 2]) + pcenter[1]

    warpR = cv2.getPerspectiveTransform(org, dst)

    dst1, dst2, dst3, dst4 = dst
    r1 = int(min(dst1[1], dst2[1]))
    r2 = int(max(dst3[1], dst4[1]))
    c1 = int(min(dst1[0], dst3[0]))
    c2 = int(max(dst2[0], dst4[0]))

    try:
        ratio = min(1.0 * h / (r2 - r1), 1.0 * w / (c2 - c1))

        dx = -c1
        dy = -r1
        T1 = np.float32([[1.0, 0, dx], [0, 1.0, dy], [0, 0, 1.0 / ratio]])
        ret = T1.dot(warpR)
    except:
        ratio = 1.0
        T1 = np.float32([[1.0, 0, 0], [0, 1.0, 0], [0, 0, 1.0]])
        ret = T1
    return ret, (-r1, -c1), ratio, dst