    def forward(
        self,
        x: torch.Tensor,
        tgt_sizes: torch.Tensor | None = None,
        attn_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        if tgt_sizes is None:
            tgt_sizes = int(math.sqrt(x.size(1)))
        if self.adaptive:
            pos_embed_arr = get_2d_sincos_pos_embed(
                self.embed_dim, tgt_sizes, version=(2, 0)
            )
            pos_embed = torch.from_numpy(pos_embed_arr).to(
                device=x.device, dtype=x.dtype
            )
        else:
            pos_embed = get_abs_pos(self.pos_embed, tgt_sizes).to(
                device=x.device, dtype=x.dtype
            )

        x, _ = self.kv_proj(x)
        x = self.ln_kv(x).permute(1, 0, 2)

        N = x.shape[1]
        q = self.ln_q(self.query)
        out = self.attn(
            self._repeat(q, N) + self.pos_embed.unsqueeze(1),
            x + pos_embed.unsqueeze(1),
            x,
            attn_mask=attn_mask,
        )[0]
        x = out.permute(1, 0, 2)
        if self.do_post_projection:
            x = self.ln_post(x)
            x = x @ self.proj
        return x