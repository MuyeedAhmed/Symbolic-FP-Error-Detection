    def demux(self, x):
        """[num_buckets, multiplex_count, ...] -> [N_obj, ...]"""
        out_shape = (self.total_valid_entries,) + x.shape[2:]
        flat = x.reshape(self.num_buckets * self.multiplex_count, -1)
        return (self.demux_matrix.to(device=x.device, dtype=x.dtype) @ flat).view(out_shape)