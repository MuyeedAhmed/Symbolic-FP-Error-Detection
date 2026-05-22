    def mux(self, x):
        """[N_obj, ...] -> [num_buckets, multiplex_count, ...]"""
        out_shape = (self.num_buckets, self.multiplex_count) + x.shape[1:]
        return (self.mux_matrix.to(device=x.device, dtype=x.dtype) @ x.reshape(self.total_valid_entries, -1)).view(out_shape)