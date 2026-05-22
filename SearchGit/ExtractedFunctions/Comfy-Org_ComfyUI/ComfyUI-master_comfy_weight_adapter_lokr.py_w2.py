    def w2(self):
        if self.w2_rebuild:
            if self.use_tucker:
                w2 = torch.einsum(
                    "i j k l, j r, i p -> p r k l",
                    self.lokr_t2,
                    self.lokr_w2_b,
                    self.lokr_w2_a,
                )
            else:
                w2 = self.lokr_w2_a @ self.lokr_w2_b
            return w2 * (self.alpha / self.rankb)
        else:
            return self.lokr_w2