    def forward(
        self,
        input_ids: torch.Tensor,
        input_values: torch.Tensor,
        attention_mask: torch.Tensor | None = None,
        padding_mask: torch.Tensor | None = None,
        return_loss: bool | None = None,
        **kwargs,
    ) -> PeAudioOutput:
        audio_outputs: BaseModelOutputWithPooling = self.audio_encoder(
            input_values=input_values, padding_mask=padding_mask, **kwargs
        )
        kwargs["output_hidden_states"] = True
        text_outputs: MaskedLMOutput = self.text_model(input_ids=input_ids, attention_mask=attention_mask, **kwargs)

        audio_embeds = audio_outputs.last_hidden_state
        audio_embeds = self.audio_head(audio_embeds)

        text_audio_embeds = text_outputs.hidden_states[-1][:, 0]
        text_audio_embeds = self.text_audio_head(text_audio_embeds)

        logits_audio_text = (audio_embeds @ text_audio_embeds.T).transpose(1, 2)
        logits_audio_text = logits_audio_text * self.text_audio_logit_scale + self.text_audio_logit_bias

        loss = None
        if return_loss:
            labels = torch.eye(logits_audio_text.shape[0], device=logits_audio_text.device)
            loss = -F.logsigmoid(labels * logits_audio_text).sum() / logits_audio_text.shape[0]

        return PeAudioOutput(
            logits_audio_text=logits_audio_text,
            text_audio_embeds=text_audio_embeds,
            audio_embeds=audio_embeds,
            text_outputs=text_outputs,
            audio_outputs=audio_outputs,
            loss=loss,
        )