def paste_instance_masks_v2(masks: np.ndarray, detected_boxes: np.ndarray,
                            image_height: int, image_width: int) -> np.ndarray:
  """Paste instance masks to generate the image segmentation (v2).

  Args:
    masks: a numpy array of shape [N, mask_height, mask_width] representing the
      instance masks w.r.t. the `detected_boxes`.
    detected_boxes: a numpy array of shape [N, 4] representing the reference
      bounding boxes.
    image_height: an integer representing the height of the image.
    image_width: an integer representing the width of the image.

  Returns:
    segms: a numpy array of shape [N, image_height, image_width] representing
      the instance masks *pasted* on the image canvas.
  """
  _, mask_height, mask_width = masks.shape

  segms = []
  for i, mask in enumerate(masks):
    box = detected_boxes[i, :]
    xmin = box[0]
    ymin = box[1]
    xmax = xmin + box[2]
    ymax = ymin + box[3]

    # Sample points of the cropped mask w.r.t. the image grid.
    # Note that these coordinates may fall beyond the image.
    # Pixel clipping will happen after warping.
    xmin_int = int(math.floor(xmin))
    xmax_int = int(math.ceil(xmax))
    ymin_int = int(math.floor(ymin))
    ymax_int = int(math.ceil(ymax))

    alpha = box[2] / (1.0 * mask_width)
    beta = box[3] / (1.0 * mask_height)
    # pylint: disable=invalid-name
    # Transformation from mask pixel indices to image coordinate.
    M_mask_to_image = np.array(
        [[alpha, 0, xmin],
         [0, beta, ymin],
         [0, 0, 1]],
        dtype=np.float32)
    # Transformation from image to cropped mask coordinate.
    M_image_to_crop = np.array(
        [[1, 0, -xmin_int],
         [0, 1, -ymin_int],
         [0, 0, 1]],
        dtype=np.float32)
    M = np.dot(M_image_to_crop, M_mask_to_image)
    # Compensate the half pixel offset that OpenCV has in the
    # warpPerspective implementation: the top-left pixel is sampled
    # at (0,0), but we want it to be at (0.5, 0.5).
    M = np.dot(
        np.dot(
            np.array([[1, 0, -0.5],
                      [0, 1, -0.5],
                      [0, 0, 1]], np.float32),
            M),
        np.array([[1, 0, 0.5],
                  [0, 1, 0.5],
                  [0, 0, 1]], np.float32))
    # pylint: enable=invalid-name
    cropped_mask = cv2.warpPerspective(
        mask.astype(np.float32), M,
        (xmax_int - xmin_int, ymax_int - ymin_int))
    cropped_mask = np.array(cropped_mask > 0.5, dtype=np.uint8)

    img_mask = np.zeros((image_height, image_width))
    x0 = max(min(xmin_int, image_width), 0)
    x1 = max(min(xmax_int, image_width), 0)
    y0 = max(min(ymin_int, image_height), 0)
    y1 = max(min(ymax_int, image_height), 0)
    img_mask[y0:y1, x0:x1] = cropped_mask[
        (y0 - ymin_int):(y1 - ymin_int),
        (x0 - xmin_int):(x1 - xmin_int)]

    segms.append(img_mask)

  segms = np.array(segms)
  return segms