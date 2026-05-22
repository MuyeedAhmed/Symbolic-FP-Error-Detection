    def forward(self, x):
        B, N, C = x.shape
        qkv = (
            self.qkv(x)
            .reshape([B, N, 3, self.num_heads, C // self.num_heads])
            .transpose([2, 0, 3, 1, 4])
        )
        q, k, v = qkv.unbind(0)  # make torchscript happy (cannot use tensor as tuple)

        attn = (q @ k.transpose([0, 1, 3, 2])) * self.scale

        attn = F.softmax(attn, axis=-1)
        attn = self.attn_drop(attn)

        x = (attn @ v).transpose([0, 2, 1, 3]).reshape([B, N, C])

        x = self.proj(x)
        x = self.proj_drop(x)
        return x