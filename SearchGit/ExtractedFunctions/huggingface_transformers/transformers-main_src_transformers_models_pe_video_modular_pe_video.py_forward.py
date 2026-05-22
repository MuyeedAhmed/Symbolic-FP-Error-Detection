    def forward(
        self,
        input_ids: torch.Tensor,
        pixel_values_videos: torch.Tensor,
        attention_mask: torch.Tensor | None = None,
        padding_mask_videos: torch.Tensor | None = None,
        return_loss: bool | None = None,
        **kwargs,
    ) -> PeVideoOutput:
        video_outputs: BaseModelOutputWithPooling = self.video_encoder(
            pixel_values_videos=pixel_values_videos, padding_mask_videos=padding_mask_videos, **kwargs
        )
        kwargs["output_hidden_states"] = True
        text_outputs: MaskedLMOutput = self.text_model(input_ids=input_ids, attention_mask=attention_mask, **kwargs)

        video_embeds = video_outputs.pooler_output
        video_embeds = self.video_head(video_embeds)

        text_video_embeds = text_outputs.hidden_states[-1][:, 0]
        text_video_embeds = self.text_video_head(text_video_embeds)

        logits_video_text = video_embeds @ text_video_embeds.T
        logits_video_text = logits_video_text * self.text_video_logit_scale + self.text_video_logit_bias

        loss = None
        if return_loss:
            labels = torch.eye(logits_video_text.shape[0], device=logits_video_text.device)
            loss = -F.logsigmoid(labels * logits_video_text).sum() / logits_video_text.shape[0]

        return PeVideoOutput(
            logits_video_text=logits_video_text,
            text_video_embeds=text_video_embeds,
            video_embeds=video_embeds,
            text_outputs=text_outputs,
            video_outputs=video_outputs,
            loss=loss,
        )