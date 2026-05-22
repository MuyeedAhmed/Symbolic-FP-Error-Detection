            def forward(self, x):
                U, S = torch.linalg.svd(x)[:2]
                reduced = U[:, :, : self.k] @ torch.diag_embed(S[:, : self.k])
                return reduced