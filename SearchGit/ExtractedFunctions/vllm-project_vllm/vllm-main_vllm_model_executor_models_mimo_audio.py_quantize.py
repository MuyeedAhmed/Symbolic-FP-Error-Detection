    def quantize(self, x):
        embed = self.embed.t()
        dist_val = -(
            x.pow(2).sum(1, keepdim=True)
            - 2 * x @ embed
            + embed.pow(2).sum(0, keepdim=True)
        )
        embed_ind = dist_val.max(dim=-1).indices
        return embed_ind