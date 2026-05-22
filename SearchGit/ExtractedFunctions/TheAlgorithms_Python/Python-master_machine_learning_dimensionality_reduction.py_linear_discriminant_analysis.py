def linear_discriminant_analysis(
    features: np.ndarray, labels: np.ndarray, classes: int, dimensions: int
) -> np.ndarray:
    """
    Linear Discriminant Analysis.

    For more details, see: https://en.wikipedia.org/wiki/Linear_discriminant_analysis.
    Parameters:
        * features: the features extracted from the dataset
        * labels: the class labels of the features
        * classes: the number of classes present in the dataset
        * dimensions: to filter the projected data for the desired dimension

    >>> test_linear_discriminant_analysis()
    """

    # Check if the dimension desired is less than the number of classes
    assert classes > dimensions

    # Check if features have been already loaded
    if features.any:
        _, eigenvectors = eigh(
            covariance_between_classes(features, labels, classes),
            covariance_within_classes(features, labels, classes),
        )
        filtered_eigenvectors = eigenvectors[:, ::-1][:, :dimensions]
        svd_matrix, _, _ = np.linalg.svd(filtered_eigenvectors)
        filtered_svd_matrix = svd_matrix[:, 0:dimensions]
        projected_data = np.dot(filtered_svd_matrix.T, features)
        logging.info("Linear Discriminant Analysis computed")

        return projected_data
    else:
        logging.basicConfig(level=logging.ERROR, format="%(message)s", force=True)
        logging.error("Dataset empty")
        raise AssertionError