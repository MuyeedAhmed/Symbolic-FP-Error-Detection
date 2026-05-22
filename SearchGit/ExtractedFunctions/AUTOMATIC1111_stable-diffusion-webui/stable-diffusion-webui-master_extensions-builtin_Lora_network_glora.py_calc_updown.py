    def calc_updown(self, orig_weight):
        w1a = self.w1a.to(orig_weight.device)
        w1b = self.w1b.to(orig_weight.device)
        w2a = self.w2a.to(orig_weight.device)
        w2b = self.w2b.to(orig_weight.device)

        output_shape = [w1a.size(0), w1b.size(1)]
        updown = ((w2b @ w1b) + ((orig_weight.to(dtype = w1a.dtype) @ w2a) @ w1a))

        return self.finalize_updown(updown, orig_weight, output_shape)