    def _build_graph(self):
        """Matrix representing a fully connected graph between each sample."""
        if self.kernel == "knn":
            self.nn_fit = None
        affinity_matrix = self._get_kernel(self.X_)
        normalizer = affinity_matrix.sum(axis=1)
        # handle spmatrix (make normalizer 1D)
        if sparse.isspmatrix(affinity_matrix):
            normalizer = np.ravel(normalizer)

        if SCIPY_VERSION_BELOW_1_12:
            if sparse.issparse(affinity_matrix):
                inv_normalizer = sparse.diags(1.0 / normalizer)
                affinity_matrix = inv_normalizer @ affinity_matrix
            else:  # Dense affinity_matrix
                affinity_matrix /= normalizer[:, np.newaxis]
            return affinity_matrix

        # same syntax for sparse or dense
        affinity_matrix /= normalizer[:, np.newaxis]
        return affinity_matrix