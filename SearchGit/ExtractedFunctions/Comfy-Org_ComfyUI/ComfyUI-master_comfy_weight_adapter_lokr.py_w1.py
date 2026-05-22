    def w1(self):
        if self.w1_rebuild:
            return (self.lokr_w1_a @ self.lokr_w1_b) * (self.alpha / self.ranka)
        else:
            return self.lokr_w1