    def forward(self, x):
        h = x.view(x.shape[0], 8, 8)
        attn = h * torch.sigmoid(h)
        # attn is both returned (user output) and needed by backward (mul saves it)
        scaled = attn * self.scale
        flat = scaled.flatten(1)
        result = flat @ self.W
        return result, attn