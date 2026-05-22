    def merge_subcluster(self, nominee_cluster, threshold):
        """Check if a cluster is worthy enough to be merged. If
        yes then merge.
        """
        new_ss = self.squared_sum_ + nominee_cluster.squared_sum_
        new_ls = self.linear_sum_ + nominee_cluster.linear_sum_
        new_n = self.n_samples_ + nominee_cluster.n_samples_
        new_centroid = (1 / new_n) * new_ls
        new_sq_norm = np.dot(new_centroid, new_centroid)

        # The squared radius of the cluster is defined:
        #   r^2  = sum_i ||x_i - c||^2 / n
        # with x_i the n points assigned to the cluster and c its centroid:
        #   c = sum_i x_i / n
        # This can be expanded to:
        #   r^2 = sum_i ||x_i||^2 / n - 2 < sum_i x_i / n, c> + n ||c||^2 / n
        # and therefore simplifies to:
        #   r^2 = sum_i ||x_i||^2 / n - ||c||^2
        sq_radius = new_ss / new_n - new_sq_norm

        if sq_radius <= threshold**2:
            (
                self.n_samples_,
                self.linear_sum_,
                self.squared_sum_,
                self.centroid_,
                self.sq_norm_,
            ) = (new_n, new_ls, new_ss, new_centroid, new_sq_norm)
            return True
        return False