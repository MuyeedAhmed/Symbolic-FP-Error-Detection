    def forward(self, x):
        B, C, H, W = x.shape
        k, v = x, x

        # calculate key vector
        features = []
        for i in range(0, len(self.k_encoder)):
            k = self.k_encoder[i](k)
            features.append(k)
        for i in range(0, len(self.k_decoder) - 1):
            k = self.k_decoder[i](k)
            # print(k.shape, features[len(self.k_decoder) - 2 - i].shape)
            k = k + features[len(self.k_decoder) - 2 - i]
        k = self.k_decoder[-1](k)

        # calculate query vector
        # TODO q=f(q,k)
        zeros = paddle.zeros((B, self.max_length, C), dtype=x.dtype)  # (B, N, C)
        q = self.pos_encoder(zeros)  # (B, N, C)
        q = self.project(q)  # (B, N, C)

        # calculate attention
        attn_scores = q @ k.flatten(2)  # (B, N, (H*W))
        attn_scores = attn_scores / (C**0.5)
        attn_scores = F.softmax(attn_scores, axis=-1)

        v = v.flatten(2).transpose([0, 2, 1])  # (B, (H*W), C)
        attn_vecs = attn_scores @ v  # (B, N, C)

        return attn_vecs, attn_scores.reshape([0, self.max_length, H, W])