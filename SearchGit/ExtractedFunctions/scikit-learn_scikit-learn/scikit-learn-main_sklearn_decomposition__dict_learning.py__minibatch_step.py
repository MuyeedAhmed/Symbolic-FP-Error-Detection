    def _minibatch_step(self, X, dictionary, random_state, step):
        """Perform the update on the dictionary for one minibatch."""
        batch_size = X.shape[0]

        # Compute code for this batch
        code = _sparse_encode(
            X,
            dictionary,
            algorithm=self._fit_algorithm,
            alpha=self.alpha,
            n_jobs=self.n_jobs,
            positive=self.positive_code,
            max_iter=self.transform_max_iter,
            verbose=self.verbose,
        )

        batch_cost = (
            0.5 * ((X - code @ dictionary) ** 2).sum()
            + self.alpha * np.sum(np.abs(code))
        ) / batch_size

        # Update inner stats
        self._update_inner_stats(X, code, batch_size, step)

        # Update dictionary
        _update_dict(
            dictionary,
            X,
            code,
            self._A,
            self._B,
            verbose=self.verbose,
            random_state=random_state,
            positive=self.positive_dict,
        )

        return batch_cost