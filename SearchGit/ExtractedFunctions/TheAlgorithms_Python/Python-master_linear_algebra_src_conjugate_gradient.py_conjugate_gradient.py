def conjugate_gradient(
    spd_matrix: np.ndarray,
    load_vector: np.ndarray,
    max_iterations: int = 1000,
    tol: float = 1e-8,
) -> Any:
    """
    Returns solution to the linear system np.dot(spd_matrix, x) = b.

    Input:
    spd_matrix is an NxN Symmetric Positive Definite (SPD) matrix.
    load_vector is an Nx1 vector.

    Output:
    x is an Nx1 vector that is the solution vector.

    >>> import numpy as np
    >>> spd_matrix = np.array([
    ... [8.73256573, -5.02034289, -2.68709226],
    ... [-5.02034289,  3.78188322,  0.91980451],
    ... [-2.68709226,  0.91980451,  1.94746467]])
    >>> b = np.array([
    ... [-5.80872761],
    ... [ 3.23807431],
    ... [ 1.95381422]])
    >>> conjugate_gradient(spd_matrix, b)
    array([[-0.63114139],
           [-0.01561498],
           [ 0.13979294]])
    """
    # Ensure proper dimensionality.
    assert np.shape(spd_matrix)[0] == np.shape(spd_matrix)[1]
    assert np.shape(load_vector)[0] == np.shape(spd_matrix)[0]
    assert _is_matrix_spd(spd_matrix)

    # Initialize solution guess, residual, search direction.
    x0 = np.zeros((np.shape(load_vector)[0], 1))
    r0 = np.copy(load_vector)
    p0 = np.copy(r0)

    # Set initial errors in solution guess and residual.
    error_residual = 1e9
    error_x_solution = 1e9
    error = 1e9

    # Set iteration counter to threshold number of iterations.
    iterations = 0

    while error > tol:
        # Save this value so we only calculate the matrix-vector product once.
        w = np.dot(spd_matrix, p0)

        # The main algorithm.

        # Update search direction magnitude.
        alpha = np.dot(r0.T, r0) / np.dot(p0.T, w)
        # Update solution guess.
        x = x0 + alpha * p0
        # Calculate new residual.
        r = r0 - alpha * w
        # Calculate new Krylov subspace scale.
        beta = np.dot(r.T, r) / np.dot(r0.T, r0)
        # Calculate new A conjuage search direction.
        p = r + beta * p0

        # Calculate errors.
        error_residual = np.linalg.norm(r - r0)
        error_x_solution = np.linalg.norm(x - x0)
        error = np.maximum(error_residual, error_x_solution)

        # Update variables.
        x0 = np.copy(x)
        r0 = np.copy(r)
        p0 = np.copy(p)

        # Update number of iterations.
        iterations += 1
        if iterations > max_iterations:
            break

    return x