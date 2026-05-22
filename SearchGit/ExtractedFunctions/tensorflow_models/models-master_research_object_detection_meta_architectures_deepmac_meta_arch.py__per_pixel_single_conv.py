def _per_pixel_single_conv(input_tensor, params, channels):
  """Convolve the given input with the given params.

  Args:
    input_tensor: A [num_instances, height, width, channels] shaped
      float tensor.
    params: A [num_instances, num_params] shaped float tensor.
    channels: int, number of channels in the convolution.

  Returns:
    output: A float tensor of shape [num_instances, height, width, channels]
  """

  input_channels = input_tensor.get_shape().as_list()[3]
  weights = params[:, :(input_channels * channels)]
  biases = params[:, (input_channels * channels):]
  num_instances = tf.shape(params)[0]

  weights = tf.reshape(weights, (num_instances, input_channels, channels))
  output = (input_tensor[:, :, tf.newaxis, :] @
            weights[:, tf.newaxis, tf.newaxis, :, :])

  output = output[:, :, 0, :, :]
  output = output + biases[:, tf.newaxis, tf.newaxis, :]
  return output