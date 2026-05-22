    def forward(
        self,
        vision_hidden_state: torch.FloatTensor,
        text_hidden_state: torch.FloatTensor,
        text_token_mask: torch.BoolTensor,
    ) -> torch.FloatTensor:
        res = vision_hidden_state @ text_hidden_state.transpose(-1, -2)
        res = res / math.sqrt(vision_hidden_state.shape[-1])
        res = res + self.bias
        res.masked_fill_(~text_token_mask[:, None, :], float("-inf"))

        # padding to max_text_len
        new_res = torch.full((*res.shape[:-1], self.max_text_len), float("-inf"), device=res.device)
        new_res[..., : res.shape[-1]] = res

        return new_res