def gaussian_filter(
    times: np.ndarray,
    values: np.ndarray,
    window_location,
) -> float:
    stdev = (max(times) - min(times)) / 2
    gaussian_distribution = scipy.stats.norm(window_location, stdev)

    coefficients = gaussian_distribution.pdf(times)
    normalized_coefficients = coefficients / sum(coefficients)
    return np.dot(values, normalized_coefficients)