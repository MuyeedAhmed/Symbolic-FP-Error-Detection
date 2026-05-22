def compute_exponential_coeffs(s: torch.Tensor, t: torch.Tensor, solver_order: int, tau_t: float) -> torch.Tensor:
    """Compute (1 + tau^2) * integral of exp((1 + tau^2) * x) * x^p dx from s to t with exp((1 + tau^2) * t) factored out, using integration by parts.

    Integral of exp((1 + tau^2) * x) * x^p dx
        = product_terms[p] - (p / (1 + tau^2)) * integral of exp((1 + tau^2) * x) * x^(p-1) dx,
    with base case p=0 where integral equals product_terms[0].

    where
        product_terms[p] = x^p * exp((1 + tau^2) * x) / (1 + tau^2).

    Construct a recursive coefficient matrix following the above recursive relation to compute all integral terms up to p = (solver_order - 1).
    Return coefficients used by the SA-Solver in data prediction mode.

    Args:
        s: Start time s.
        t: End time t.
        solver_order: Current order of the solver.
        tau_t: Stochastic strength parameter in the SDE.

    Returns:
        Exponential coefficients used in data prediction, with exp((1 + tau^2) * t) factored out, ordered from p=0 to p=solver_order−1, shape (solver_order,).
    """
    tau_mul = 1 + tau_t ** 2
    h = t - s
    p = torch.arange(solver_order, dtype=s.dtype, device=s.device)

    # product_terms after factoring out exp((1 + tau^2) * t)
    # Includes (1 + tau^2) factor from outside the integral
    product_terms_factored = (t ** p - s ** p * (-tau_mul * h).exp())

    # Lower triangular recursive coefficient matrix
    # Accumulates recursive coefficients based on p / (1 + tau^2)
    recursive_depth_mat = p.unsqueeze(1) - p.unsqueeze(0)
    log_factorial = (p + 1).lgamma()
    recursive_coeff_mat = log_factorial.unsqueeze(1) - log_factorial.unsqueeze(0)
    if tau_t > 0:
        recursive_coeff_mat = recursive_coeff_mat - (recursive_depth_mat * math.log(tau_mul))
    signs = torch.where(recursive_depth_mat % 2 == 0, 1.0, -1.0)
    recursive_coeff_mat = (recursive_coeff_mat.exp() * signs).tril()

    return recursive_coeff_mat @ product_terms_factored