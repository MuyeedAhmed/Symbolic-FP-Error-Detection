def principal_component_analysis(features: np.ndarray, dimensions: int) -> np.ndarray:
    """
    Principal Component Analysis.

    For more details, see: https://en.wikipedia.org/wiki/Principal_component_analysis.
    Parameters:
        * features: the features extracted from the dataset
        * dimensions: to filter the projected data for the desired dimension

    >>> test_principal_component_analysis()
    """

    # Check if the features have been loaded
    if features.any():
        data_mean = features.mean(1)
        # Center the dataset
        centered_data = features - np.reshape(data_mean, (data_mean.size, 1))
        covariance_matrix = np.dot(centered_data, centered_data.T) / features.shape[1]
        _, eigenvectors = np.linalg.eigh(covariance_matrix)
        # Take all the columns in the reverse order (-1), and then takes only the first
        filtered_eigenvectors = eigenvectors[:, ::-1][:, 0:dimensions]
        # Project the database on the new space
        projected_data = np.dot(filtered_eigenvectors.T, features)
        logging.info("Principal Component Analysis computed")

        return projected_data
    else:
        logging.basicConfig(level=logging.ERROR, format="%(message)s", force=True)
        logging.error("Dataset empty")
        raise AssertionError