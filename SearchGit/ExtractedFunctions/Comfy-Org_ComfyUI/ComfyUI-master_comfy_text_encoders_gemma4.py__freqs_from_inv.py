    def _freqs_from_inv(self, inv_freq, position_ids, device, dtype):
        """Compute cos/sin from stored inv_freq"""
        inv_exp = inv_freq[None, :, None].float().expand(position_ids.shape[0], -1, 1).to(device)
        pos_exp = position_ids[:, None, :].float()
        freqs = (inv_exp @ pos_exp).transpose(1, 2)
        emb = torch.cat((freqs, freqs), dim=-1)
        return emb.cos().unsqueeze(1).to(dtype), emb.sin().unsqueeze(1).to(dtype)