    def _fine_matching(
        self,
        fine_features_0: torch.Tensor,
        fine_features_1: torch.Tensor,
        coarse_matched_keypoints: torch.Tensor,
        fine_scale: float,
    ) -> torch.Tensor:
        """
        For each coarse pixel with a corresponding window of fine features, compute the matching confidence between fine
        features in the first image and the second image.

        Fine features are sliced in two part :
        - The first part used for the first stage are the first fine_hidden_size - config.fine_matching_slicedim (64 - 8
         = 56 by default) features.
        - The second part used for the second stage are the last config.fine_matching_slicedim (8 by default) features.

        Each part is used to compute a fine confidence tensor of the following shape :
        (batch_size, (coarse_height * coarse_width), fine_window_size, fine_window_size)
        They correspond to the score between each fine pixel in the first image and each fine pixel in the second image.

        Args:
            fine_features_0 (`torch.Tensor` of shape `(num_matches, fine_kernel_size ** 2, fine_kernel_size ** 2)`):
                Fine features from the first image
            fine_features_1 (`torch.Tensor` of shape `(num_matches, (fine_kernel_size + 2) ** 2, (fine_kernel_size + 2)
            ** 2)`):
                Fine features from the second image
            coarse_matched_keypoints (`torch.Tensor` of shape `(2, num_matches, 2)`):
                Keypoint coordinates found in coarse matching for the first and second image
            fine_scale (`int`):
                Scale between the size of fine features and coarse features

        Returns:
            fine_coordinates (`torch.Tensor` of shape `(2, num_matches, 2)`):
                Matched keypoint between the first and the second image. All matched keypoints are concatenated in the
                second dimension.

        """
        batch_size, num_keypoints, fine_window_size, fine_embed_dim = fine_features_0.shape
        fine_matching_slice_dim = self.config.fine_matching_slice_dim

        fine_kernel_size = torch_int(fine_window_size**0.5)

        # Split fine features into first and second stage features
        split_fine_features_0 = torch.split(fine_features_0, fine_embed_dim - fine_matching_slice_dim, -1)
        split_fine_features_1 = torch.split(fine_features_1, fine_embed_dim - fine_matching_slice_dim, -1)

        # Retrieve first stage fine features
        fine_features_0 = split_fine_features_0[0]
        fine_features_1 = split_fine_features_1[0]

        # Normalize first stage fine features
        fine_features_0 = fine_features_0 / fine_features_0.shape[-1] ** 0.5
        fine_features_1 = fine_features_1 / fine_features_1.shape[-1] ** 0.5

        # Compute first stage confidence
        fine_confidence = fine_features_0 @ fine_features_1.transpose(-1, -2)
        fine_confidence = nn.functional.softmax(fine_confidence, 1) * nn.functional.softmax(fine_confidence, 2)
        fine_confidence = fine_confidence.reshape(
            batch_size, num_keypoints, fine_window_size, fine_kernel_size + 2, fine_kernel_size + 2
        )
        fine_confidence = fine_confidence[..., 1:-1, 1:-1]
        first_stage_fine_confidence = fine_confidence.reshape(
            batch_size, num_keypoints, fine_window_size, fine_window_size
        )

        fine_indices, fine_matches = self._get_first_stage_fine_matching(
            first_stage_fine_confidence,
            coarse_matched_keypoints,
            fine_window_size,
            fine_scale,
        )

        # Retrieve second stage fine features
        fine_features_0 = split_fine_features_0[1]
        fine_features_1 = split_fine_features_1[1]

        # Normalize second stage fine features
        fine_features_1 = fine_features_1 / fine_matching_slice_dim**0.5

        # Compute second stage fine confidence
        second_stage_fine_confidence = fine_features_0 @ fine_features_1.transpose(-1, -2)

        fine_coordinates = self._get_second_stage_fine_matching(
            fine_indices,
            fine_matches,
            second_stage_fine_confidence,
            fine_window_size,
            fine_scale,
        )

        return fine_coordinates