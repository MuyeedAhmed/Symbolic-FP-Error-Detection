    def forward(self, x):
        shape, dtype = x.shape, x.dtype
        x = self.preprocess(x)

        self.init_embed_(x)

        embed_ind = self.quantize(x)
        embed_onehot = F.one_hot(embed_ind, self.codebook_size).type(dtype)
        embed_ind = self.postprocess_emb(embed_ind, shape)
        quantize = self.dequantize(embed_ind)

        if self.training:
            ### Update codebook by EMA
            embed_onehot_sum = embed_onehot.sum(0)  # [cb-size,]
            embed_sum = x.t() @ embed_onehot  # [D, cb-size]
            if is_distributed():
                dist.all_reduce(embed_onehot_sum)
                dist.all_reduce(embed_sum)
            # Update ema cluster count N_i^t, eq. (6) in vqvae paper
            self.cluster_size.data.mul_(self.decay).add_(embed_onehot_sum, alpha=1 - self.decay)
            # Update ema embed: eq. (7) in vqvae paper
            self.embed_avg.data.mul_(self.decay).add_(embed_sum.t(), alpha=1 - self.decay)
            # apply laplace smoothing
            n = self.cluster_size.sum()
            cluster_size = (self.cluster_size + self.epsilon) / (n + self.codebook_size * self.epsilon) * n
            # Update ema embed: eq. (8) in vqvae paper
            embed_normalized = self.embed_avg / cluster_size.unsqueeze(1)
            self.embed.data.copy_(embed_normalized)

            # We do the expiry of code at that point as buffers are in sync
            # and all the workers will take the same decision.
            self.expire_codes_(x)

        return quantize, embed_ind