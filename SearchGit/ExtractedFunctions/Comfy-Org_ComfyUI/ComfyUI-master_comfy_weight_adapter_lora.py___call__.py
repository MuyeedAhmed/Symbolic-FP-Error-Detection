    def __call__(self, w):
        org_dtype = w.dtype
        if self.lora_mid is None:
            diff = self.lora_up.weight @ self.lora_down.weight
        else:
            diff = tucker_weight_from_conv(
                self.lora_up.weight, self.lora_down.weight, self.lora_mid.weight
            )
        scale = self.alpha / self.rank
        weight = w + scale * diff.reshape(w.shape)
        return weight.to(org_dtype)