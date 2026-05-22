    def forward(self, x):

        B, H, W, _ = x.shape
        qkv = (
            self.qkv(x)
            .reshape([B, H * W, 3, self.num_heads, -1])
            .transpose([2, 0, 3, 1, 4])
        )
        q, k, v = qkv.reshape([3, B * self.num_heads, H * W, -1]).unbind(0)
        attn = (q * self.scale) @ k.transpose([0, 2, 1])

        if self.use_rel_pos:
            attn = add_decomposed_rel_pos(
                attn, q, self.rel_pos_h, self.rel_pos_w, (H, W), (H, W)
            )
        attn = F.softmax(attn, axis=-1)
        x = (
            (attn @ v)
            .reshape([B, self.num_heads, H, W, -1])
            .transpose([0, 2, 3, 1, 4])
            .reshape([B, H, W, -1])
        )
        x = self.proj(x)

        return x