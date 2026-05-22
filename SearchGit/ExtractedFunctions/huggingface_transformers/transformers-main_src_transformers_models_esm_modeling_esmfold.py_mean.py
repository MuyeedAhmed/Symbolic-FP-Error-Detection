    def mean(self):
        return (self.logits.softmax(-1) @ self.v_bins.unsqueeze(1)).squeeze(-1)