    def forward(
        self,
        vision_hidden_state: torch.FloatTensor,
        text_hidden_state: torch.FloatTensor,
        text_token_mask: torch.BoolTensor,
    ) -> torch.FloatTensor:
        output = vision_hidden_state @ text_hidden_state.transpose(-1, -2)
        output = output.masked_fill(~text_token_mask[:, None, :], float("-inf"))

        # padding to max_text_len
        new_output = torch.full((*output.shape[:-1], self.max_text_len), float("-inf"), device=output.device)
        new_output[..., : output.shape[-1]] = output

        return new_output