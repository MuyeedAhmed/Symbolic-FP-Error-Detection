    def forward(
        self,
        image_embeddings: torch.Tensor,
        image_positional_embeddings: torch.Tensor,
        sparse_prompt_embeddings: torch.Tensor,
        dense_prompt_embeddings: torch.Tensor,
        multimask_output: bool,
        high_resolution_features: list[torch.Tensor],
        attention_similarity: torch.Tensor | None = None,
        target_embedding: torch.Tensor | None = None,
        **kwargs: Unpack[TransformersKwargs],
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        Predict masks given image and prompt embeddings.

        Args:
            image_embeddings (`torch.Tensor`):
                The embeddings from the image encoder.
            image_positional_embeddings (`torch.Tensor`):
                Positional encoding with the shape of image_embeddings.
            sparse_prompt_embeddings (`torch.Tensor`):
                The embeddings of the points and boxes.
            dense_prompt_embeddings (`torch.Tensor`):
                The embeddings of the mask inputs.
            multimask_output (`bool`):
                Whether to return multiple masks or a single mask.
            high_resolution_features (`list[torch.Tensor]`, *optional*):
                The high-resolution features from the vision encoder.
            attention_similarity (`torch.Tensor`, *optional*):
                The attention similarity tensor.
            target_embedding (`torch.Tensor`, *optional*):
                The target embedding.
        """
        batch_size, num_channels, height, width = image_embeddings.shape
        point_batch_size = sparse_prompt_embeddings.shape[1]
        # Concatenate output tokens
        output_tokens = torch.cat(
            [
                self.obj_score_token.weight,
                self.iou_token.weight,
                self.mask_tokens.weight,
            ],
            dim=0,
        )
        output_tokens = output_tokens.repeat(batch_size, point_batch_size, 1, 1)

        if sparse_prompt_embeddings.shape[0] != 0:
            tokens = torch.cat((output_tokens, sparse_prompt_embeddings), dim=2)
        else:
            tokens = output_tokens
        point_embeddings = tokens.to(self.iou_token.weight.dtype)

        # Expand per-image data in batch direction to be per-mask
        image_embeddings = image_embeddings + dense_prompt_embeddings
        image_embeddings = image_embeddings.repeat_interleave(point_batch_size, dim=0)
        image_positional_embeddings = image_positional_embeddings.repeat_interleave(point_batch_size, 0)
        # Run the transformer
        point_embeddings, image_embeddings = self.transformer(
            point_embeddings=point_embeddings,
            image_embeddings=image_embeddings,
            image_positional_embeddings=image_positional_embeddings,
            attention_similarity=attention_similarity,
            target_embedding=target_embedding,
            **kwargs,
        )
        iou_token_out = point_embeddings[:, :, 1, :]
        mask_tokens_out = point_embeddings[:, :, 2 : (2 + self.num_mask_tokens), :]

        # Upscale mask embeddings and predict masks using the mask tokens
        image_embeddings = image_embeddings.transpose(2, 3).view(
            batch_size * point_batch_size, num_channels, height, width
        )

        feat_s0, feat_s1 = high_resolution_features
        feat_s0 = feat_s0.repeat_interleave(point_batch_size, dim=0)
        feat_s1 = feat_s1.repeat_interleave(point_batch_size, dim=0)
        upscaled_embedding = self.upscale_conv1(image_embeddings) + feat_s1
        upscaled_embedding = self.activation(self.upscale_layer_norm(upscaled_embedding))
        upscaled_embedding = self.activation(self.upscale_conv2(upscaled_embedding) + feat_s0)

        hyper_in_list: list[torch.Tensor] = []
        for i in range(self.num_mask_tokens):
            current_mlp = self.output_hypernetworks_mlps[i]
            hyper_in_list += [current_mlp(mask_tokens_out[:, :, i, :])]
        hyper_in = torch.stack(hyper_in_list, dim=2)

        _, num_channels, height, width = upscaled_embedding.shape
        upscaled_embedding = upscaled_embedding.view(batch_size, point_batch_size, num_channels, height * width)
        masks = (hyper_in @ upscaled_embedding).view(batch_size, point_batch_size, -1, height, width)

        # Generate mask quality predictions
        iou_pred = self.iou_prediction_head(iou_token_out)
        object_score_logits = self.pred_obj_score_head(point_embeddings[:, :, 0, :])

        # Select the correct mask or masks for output
        if multimask_output:
            mask_slice = slice(1, None)
            masks = masks[:, :, mask_slice, :, :]
            iou_pred = iou_pred[:, :, mask_slice]
        elif self.dynamic_multimask_via_stability and not self.training:
            mask_slice = slice(0, 1)
            masks, iou_pred = self._dynamic_multimask_via_stability(masks, iou_pred)
        else:
            mask_slice = slice(0, 1)
            masks = masks[:, :, mask_slice, :, :]
            iou_pred = iou_pred[:, :, mask_slice]

        sam_tokens_out = mask_tokens_out[:, :, mask_slice]  # [b, 3, c] shape

        return masks, iou_pred, sam_tokens_out, object_score_logits