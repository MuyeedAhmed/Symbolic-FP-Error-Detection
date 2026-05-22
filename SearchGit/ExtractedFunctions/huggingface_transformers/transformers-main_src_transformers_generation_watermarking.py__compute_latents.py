    def _compute_latents(self, g_values: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor]:
        """Computes the unique token probability distribution given g-values.

        Args:
            g_values (`torch.Tensor` of shape `(batch_size, seq_len, watermarking_depth)`):
                PRF values.

        Returns:
            p_one_unique_token and p_two_unique_tokens, both of shape
            [batch_size, seq_len, watermarking_depth]. p_one_unique_token[i,t,l]
            gives the probability of there being one unique token in a tournament
            match on layer l, on timestep t, for batch item i.
            p_one_unique_token[i,t,l] + p_two_unique_token[i,t,l] = 1.
        """
        # Tile g-values to produce feature vectors for predicting the latents
        # for each layer in the tournament; our model for the latents psi is a
        # logistic regression model psi = sigmoid(delta * x + beta).

        # [batch_size, seq_len, watermarking_depth, watermarking_depth]
        x = torch.repeat_interleave(torch.unsqueeze(g_values, dim=-2), self.watermarking_depth, axis=-2)

        # mask all elements above -1 diagonal for autoregressive factorization
        x = torch.tril(x, diagonal=-1)

        # [batch_size, seq_len, watermarking_depth]
        # (i, j, k, l) x (i, j, k, l) -> (i, j, k) einsum equivalent
        logits = (self.delta[..., None, :] @ x.type(self.delta.dtype)[..., None]).squeeze() + self.beta

        p_two_unique_tokens = torch.sigmoid(logits)
        p_one_unique_token = 1 - p_two_unique_tokens
        return p_one_unique_token, p_two_unique_tokens