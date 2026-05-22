    def forward(self, descriptors: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        batch_size, num_keypoints, descriptor_dim = descriptors.shape
        # Final projection and similarity computation
        m_descriptors = self.final_projection(descriptors)
        m_descriptors = m_descriptors / torch.tensor(self.descriptor_dim, device=m_descriptors.device) ** 0.25
        m_descriptors = m_descriptors.reshape(batch_size // 2, 2, num_keypoints, descriptor_dim)
        m_descriptors0 = m_descriptors[:, 0]
        m_descriptors1 = m_descriptors[:, 1]
        similarity = m_descriptors0 @ m_descriptors1.transpose(-1, -2)
        if mask is not None:
            mask = mask.reshape(batch_size // 2, 2, num_keypoints)
            mask0 = mask[:, 0].unsqueeze(-1)
            mask1 = mask[:, 1].unsqueeze(-1).transpose(-1, -2)
            mask = mask0 * mask1
            similarity = similarity.masked_fill(mask == 0, torch.finfo(similarity.dtype).min)

        # Compute matchability of descriptors
        matchability = self.matchability(descriptors)
        matchability = matchability.reshape(batch_size // 2, 2, num_keypoints, 1)
        matchability_0 = matchability[:, 0]
        matchability_1 = matchability[:, 1]

        # Compute scores from similarity and matchability
        scores = sigmoid_log_double_softmax(similarity, matchability_0, matchability_1)
        return scores