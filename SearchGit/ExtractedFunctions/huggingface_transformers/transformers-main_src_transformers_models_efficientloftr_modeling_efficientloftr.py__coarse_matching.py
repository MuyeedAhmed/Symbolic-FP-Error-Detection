    def _coarse_matching(
        self, coarse_features: torch.Tensor, coarse_scale: float
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        """
        For each image pair, compute the matching confidence between each coarse element (by default (image_height / 8)
        * (image_width / 8 elements)) from the first image to the second image.

        Note:
            This step can be done as a postprocessing step, because does not involve any model weights/params.
            However, we keep it in the modeling code for consistency with other keypoint matching models AND for
            easier torch.compile/torch.export (all ops are in torch).

        Args:
            coarse_features (`torch.Tensor` of shape `(batch_size, 2, hidden_size, coarse_height, coarse_width)`):
                Coarse features
            coarse_scale (`float`): Scale between the image size and the coarse size

        Returns:
            keypoints (`torch.Tensor` of shape `(batch_size, 2, num_matches, 2)`):
                Keypoints coordinates.
            matching_scores (`torch.Tensor` of shape `(batch_size, 2, num_matches)`):
                The confidence matching score of each keypoint.
            matched_indices (`torch.Tensor` of shape `(batch_size, 2, num_matches)`):
                Indices which indicates which keypoint in an image matched with which keypoint in the other image. For
                both image in the pair.
        """
        batch_size, _, embed_dim, height, width = coarse_features.shape

        # (batch_size, 2, embed_dim, height, width) -> (batch_size, 2, height * width, embed_dim)
        coarse_features = coarse_features.permute(0, 1, 3, 4, 2)
        coarse_features = coarse_features.reshape(batch_size, 2, -1, embed_dim)

        coarse_features = coarse_features / coarse_features.shape[-1] ** 0.5
        coarse_features_0 = coarse_features[:, 0]
        coarse_features_1 = coarse_features[:, 1]

        similarity = coarse_features_0 @ coarse_features_1.transpose(-1, -2)
        similarity = similarity / self.config.coarse_matching_temperature

        if self.config.coarse_matching_skip_softmax:
            confidence = similarity
        else:
            confidence = nn.functional.softmax(similarity, 1) * nn.functional.softmax(similarity, 2)

        confidence = confidence.view(batch_size, height, width, height, width)
        matched_indices, matching_scores = self._get_matches_from_scores(confidence)

        keypoints = torch.stack([matched_indices % width, matched_indices // width], dim=-1) * coarse_scale

        return keypoints, matching_scores, matched_indices