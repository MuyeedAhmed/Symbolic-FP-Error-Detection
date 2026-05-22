def rotate_bbox(img, text_polys, angle, scale=1):
    """
    from https://github.com/WenmuZhou/DBNet.pytorch/blob/master/data_loader/modules/augment.py
    Args:
        img: np.ndarray
        text_polys: np.ndarray N*4*2
        angle: int
        scale: int

    Returns:

    """
    w = img.shape[1]
    h = img.shape[0]

    rangle = np.deg2rad(angle)
    nw = abs(np.sin(rangle) * h) + abs(np.cos(rangle) * w)
    nh = abs(np.cos(rangle) * h) + abs(np.sin(rangle) * w)
    rot_mat = cv2.getRotationMatrix2D((nw * 0.5, nh * 0.5), angle, scale)
    rot_move = np.dot(rot_mat, np.array([(nw - w) * 0.5, (nh - h) * 0.5, 0]))
    rot_mat[0, 2] += rot_move[0]
    rot_mat[1, 2] += rot_move[1]

    # ---------------------- rotate box ----------------------
    rot_text_polys = list()
    for bbox in text_polys:
        point1 = np.dot(rot_mat, np.array([bbox[0, 0], bbox[0, 1], 1]))
        point2 = np.dot(rot_mat, np.array([bbox[1, 0], bbox[1, 1], 1]))
        point3 = np.dot(rot_mat, np.array([bbox[2, 0], bbox[2, 1], 1]))
        point4 = np.dot(rot_mat, np.array([bbox[3, 0], bbox[3, 1], 1]))
        rot_text_polys.append([point1, point2, point3, point4])
    return np.array(rot_text_polys, dtype=np.float32)