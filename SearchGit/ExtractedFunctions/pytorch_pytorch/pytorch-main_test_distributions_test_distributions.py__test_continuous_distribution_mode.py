    def _test_continuous_distribution_mode(self, dist, sanitized_mode, batch_isfinite):
        # We perturb the mode in the unconstrained space and expect the log probability to decrease.
        num_points = 10
        transform = transform_to(dist.support)
        unconstrained_mode = transform.inv(sanitized_mode)
        perturbation = 1e-5 * (
            torch.rand((num_points,) + unconstrained_mode.shape) - 0.5
        )
        perturbed_mode = transform(perturbation + unconstrained_mode)
        log_prob_mode = dist.log_prob(sanitized_mode)
        log_prob_other = dist.log_prob(perturbed_mode)
        delta = log_prob_mode - log_prob_other

        # We pass the test with a small tolerance to allow for rounding and manually set the
        # difference to zero if both log probs are infinite with the same sign.
        both_infinite_with_same_sign = (log_prob_mode == log_prob_other) & (
            log_prob_mode.abs() == inf
        )
        delta[both_infinite_with_same_sign] = 0.0
        ordering = (delta > -1e-12).all(axis=0)
        self.assertTrue(ordering[batch_isfinite].all())