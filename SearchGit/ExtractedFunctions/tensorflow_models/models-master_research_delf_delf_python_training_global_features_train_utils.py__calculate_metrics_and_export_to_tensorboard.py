def _calculate_metrics_and_export_to_tensorboard(vecs, qvecs, dataset, cfg,
                                                 writer, epoch, whiten=False):
  """
  Calculates metrics and exports them to tensorboard.

  Args:
    vecs: Numpy array dataset global descriptors.
    qvecs: Numpy array query global descriptors.
    dataset: String, one of `_TEST_DATASET_NAMES`.
    cfg: Dataset configuration.
    writer: Tensorboard writer.
    epoch: Integer, epoch number.
    whiten: Boolean, whether the metrics are with for whitening used as a
      post-processing step. Affects the name of the extracted TensorBoard
      metrics.
  """
  # Search, rank and print test set metrics.
  scores = np.dot(vecs.T, qvecs)
  ranks = np.transpose(np.argsort(-scores, axis=0))

  metrics = global_features_utils.compute_metrics_and_print(dataset, ranks,
                                                            cfg['gnd'])
  # Save calculated metrics in a tensorboard format.
  if writer:
    if whiten:
      metric_names = ['test_accuracy_whiten_{}_E'.format(dataset),
                      'test_accuracy_whiten_{}_M'.format(dataset),
                      'test_accuracy_whiten_{}_H'.format(dataset)]
    else:
      metric_names = ['test_accuracy_{}_E'.format(dataset),
                      'test_accuracy_{}_M'.format(dataset),
                      'test_accuracy_{}_H'.format(dataset)]
    tf.summary.scalar(metric_names[0], metrics[0][0], step=epoch)
    tf.summary.scalar(metric_names[1], metrics[1][0], step=epoch)
    tf.summary.scalar(metric_names[2], metrics[2][0], step=epoch)
    writer.flush()
  return None