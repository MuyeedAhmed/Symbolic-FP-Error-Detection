    def forward(self, image_embeddings, image_pe, sparse_prompt_embeddings, dense_prompt_embeddings,
                high_res_features=None, multimask_output=False, return_all=False):
        B = sparse_prompt_embeddings.shape[0]
        ref = sparse_prompt_embeddings
        # Token order: [obj_score(1), iou(1), mask(num_mask_tokens)]
        tokens = torch.cat([cast_to_input(self.obj_score_token.weight, ref),
                            cast_to_input(self.iou_token.weight, ref),
                            cast_to_input(self.mask_tokens.weight, ref)], dim=0)
        tokens = torch.cat([tokens.unsqueeze(0).expand(B, -1, -1), sparse_prompt_embeddings], dim=1)

        src = image_embeddings
        if src.shape[0] != B:
            src = src.expand(B, -1, -1, -1)
        src = src + dense_prompt_embeddings
        pos_src = image_pe.expand(B, -1, -1, -1)

        b, c, h, w = src.shape
        src_flat = src.flatten(2).permute(0, 2, 1)
        pos_flat = pos_src.flatten(2).permute(0, 2, 1)

        hs, src_out = self.transformer(src_flat, pos_flat, tokens)

        obj_score_token_out = hs[:, 0, :]
        iou_token_out = hs[:, 1, :]
        mask_tokens_out = hs[:, 2:2 + self.num_mask_tokens, :]

        src_out = src_out.permute(0, 2, 1).view(b, c, h, w)
        upscaled = _upscale_masks(self.output_upscaling, self.conv_s0, self.conv_s1, src_out, high_res_features)

        hyper_in = torch.stack([
            mlp(mask_tokens_out[:, i, :]) for i, mlp in enumerate(self.output_hypernetworks_mlps)
        ], dim=1)

        masks = (hyper_in @ upscaled.flatten(2)).view(B, self.num_mask_tokens, upscaled.shape[2], upscaled.shape[3])
        iou_pred = self.iou_prediction_head(iou_token_out)
        object_score_logits = self.pred_obj_score_head(obj_score_token_out)

        if multimask_output:
            out_masks = masks[:, 1:]
            out_iou = iou_pred[:, 1:]
            out_tokens = mask_tokens_out[:, 1:]
        else:
            out_masks = masks[:, 0:1]
            out_iou = iou_pred[:, 0:1]
            out_tokens = mask_tokens_out[:, 0:1]

        if return_all:
            return out_masks, out_iou, out_tokens, object_score_logits
        return out_masks, out_iou