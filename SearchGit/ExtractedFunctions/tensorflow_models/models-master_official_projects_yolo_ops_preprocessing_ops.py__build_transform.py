def _build_transform(image,
                     perspective=0.00,
                     degrees=0.0,
                     scale_min=1.0,
                     scale_max=1.0,
                     translate=0.0,
                     random_pad=False,
                     desired_size=None,
                     seed=None):
  """Builds a unified affine transformation to spatially augment the image."""

  height, width = get_image_shape(image)
  ch = height = tf.cast(height, tf.float32)
  cw = width = tf.cast(width, tf.float32)
  deg_to_rad = lambda x: tf.cast(x, tf.float32) * np.pi / 180.0

  if desired_size is not None:
    desired_size = tf.cast(desired_size, tf.float32)
    ch = desired_size[0]
    cw = desired_size[1]

  # Compute the center of the image in the output resulution.
  center = tf.eye(3, dtype=tf.float32)
  center = tf.tensor_scatter_nd_update(center, [[0, 2], [1, 2]],
                                       [-cw / 2, -ch / 2])
  center_boxes = tf.tensor_scatter_nd_update(center, [[0, 2], [1, 2]],
                                             [cw / 2, ch / 2])

  # Compute a random rotation to apply.
  rotation = tf.eye(3, dtype=tf.float32)
  a = deg_to_rad(random_uniform_strong(-degrees, degrees, seed=seed))
  cos = tf.math.cos(a)
  sin = tf.math.sin(a)
  rotation = tf.tensor_scatter_nd_update(rotation,
                                         [[0, 0], [0, 1], [1, 0], [1, 1]],
                                         [cos, -sin, sin, cos])
  rotation_boxes = tf.tensor_scatter_nd_update(rotation,
                                               [[0, 0], [0, 1], [1, 0], [1, 1]],
                                               [cos, sin, -sin, cos])

  # Compute a random prespective change to apply.
  prespective_warp = tf.eye(3)
  px = random_uniform_strong(-perspective, perspective, seed=seed)
  py = random_uniform_strong(-perspective, perspective, seed=seed)
  prespective_warp = tf.tensor_scatter_nd_update(prespective_warp,
                                                 [[2, 0], [2, 1]], [px, py])
  prespective_warp_boxes = tf.tensor_scatter_nd_update(prespective_warp,
                                                       [[2, 0], [2, 1]],
                                                       [-px, -py])

  # Compute a random scaling to apply.
  scale = tf.eye(3, dtype=tf.float32)
  s = random_uniform_strong(scale_min, scale_max, seed=seed)
  scale = tf.tensor_scatter_nd_update(scale, [[0, 0], [1, 1]], [1 / s, 1 / s])
  scale_boxes = tf.tensor_scatter_nd_update(scale, [[0, 0], [1, 1]], [s, s])

  # Compute a random Translation to apply.
  translation = tf.eye(3)
  if (random_pad and height * s < ch and width * s < cw):
    # The image is contained within the image and arbitrarily translated to
    # locations with in the image.
    center = center_boxes = tf.eye(3, dtype=tf.float32)
    tx = random_uniform_strong(-1, 0, seed=seed) * (cw / s - width)
    ty = random_uniform_strong(-1, 0, seed=seed) * (ch / s - height)
  else:
    # The image can be translated outside of the output resolution window
    # but the image is translated relative to the output resolution not the
    # input image resolution.
    tx = random_uniform_strong(0.5 - translate, 0.5 + translate, seed=seed)
    ty = random_uniform_strong(0.5 - translate, 0.5 + translate, seed=seed)

    # Center and Scale the image such that the window of translation is
    # contained to the output resolution.
    dx, dy = (width - cw / s) / width, (height - ch / s) / height
    sx, sy = 1 - dx, 1 - dy
    bx, by = dx / 2, dy / 2
    tx, ty = bx + (sx * tx), by + (sy * ty)

    # Scale the translation to width and height of the image.
    tx *= width
    ty *= height

  translation = tf.tensor_scatter_nd_update(translation, [[0, 2], [1, 2]],
                                            [tx, ty])
  translation_boxes = tf.tensor_scatter_nd_update(translation, [[0, 2], [1, 2]],
                                                  [-tx, -ty])

  # Use repeated matric multiplications to combine all the image transforamtions
  # into a single unified augmentation operation M is applied to the image
  # Mb is to apply to the boxes. The order of matrix multiplication is
  # important. First, Translate, then Scale, then Rotate, then Center, then
  # finally alter the Prepsective.
  affine = (translation @ scale @ rotation @ center @ prespective_warp)
  affine_boxes = (
      prespective_warp_boxes @ center_boxes @ rotation_boxes @ scale_boxes
      @ translation_boxes)
  return affine, affine_boxes, s