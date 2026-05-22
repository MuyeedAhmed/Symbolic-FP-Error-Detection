    def forward(self, x):
        B, N, C = x.shape
        C = C // 3
        qkv = x.reshape([B, N, 3, C // self.head_dim, self.head_dim]).transpose(
            [2, 0, 3, 1, 4]
        )
        q, k, v = qkv.unbind(0)
        attn = q @ k.transpose([0, 1, 3, 2]) * self.scale
        attn = F.softmax(attn, -1)
        attn = self.attn_drop(attn)
        x = attn @ v
        x = x.transpose([0, 2, 1]).reshape([B, N, C])
        return x