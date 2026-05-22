    def __call__(self, w):
        org_dtype = w.dtype
        I = torch.eye(self.block_size, device=self.oft_blocks.device)

        ## generate r
        # for Q = -Q^T
        q = self.oft_blocks - self.oft_blocks.transpose(1, 2)
        normed_q = q
        if self.constraint:
            q_norm = torch.norm(q) + 1e-8
            if q_norm > self.constraint:
                normed_q = q * self.constraint / q_norm
        # use float() to prevent unsupported type
        r = (I + normed_q) @ (I - normed_q).float().inverse()

        ## Apply chunked matmul on weight
        _, *shape = w.shape
        org_weight = w.to(dtype=r.dtype)
        org_weight = org_weight.unflatten(0, (self.block_num, self.block_size))
        # Init R=0, so add I on it to ensure the output of step0 is original model output
        weight = torch.einsum(
            "k n m, k n ... -> k m ...",
            r,
            org_weight,
        ).flatten(0, 1)
        if self.rescaled:
            weight = self.rescale * weight
        return weight.to(org_dtype)