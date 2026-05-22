    def _torch_extract_fbank_features(
        self, waveform: "torch.FloatTensor", audio_lengths: "torch.Tensor", device: str = "cpu"
    ) -> "torch.FloatTensor":
        """
        Compute the log mel-scaled spectrogram of batched waveforms using PyTorch's FFT implementation.

        Args:
            waveform (torch.FloatTensor` of shape `(batch_size, max_audio_length)`):
                The batched waveforms.
            audio_lengths (`torch.Tensor` of shape `(batch_size,)`):
                The lengths of the waveforms along the max_audio_length dimension.
            device (`str`, *optional*, defaults to "cpu"):
                The device to run the computation on. (e.g., "cpu", "cuda")

        Returns:
            `torch.FloatTensor` of shape `(batch_size, max_feature_length, feature_size)`:
                The log mel-scaled spectrogram of the batched waveforms.
        """
        fft_window = torch.hamming_window(self.win_length, periodic=False, device=device, dtype=torch.float64)

        # batched implementation
        batch_size = waveform.shape[0]
        frames = waveform.unfold(-1, self.win_length, self.hop_length)

        # ---
        # the unbatched (and unpaded) original implementation skips last few audio values that can't be included in a frame
        # we need to ensure that the corresponding frames for the padded input also mask these values
        if batch_size > 1:
            frames = frames.clone()
            # concerned batch indices
            to_mask_batch_idxs = torch.arange(batch_size)[audio_lengths != audio_lengths.max()]
            if to_mask_batch_idxs.numel() > 0:
                batch_idxs_down = (audio_lengths[to_mask_batch_idxs] - self.win_length) // self.hop_length + 1
                batch_idxs_up = (audio_lengths[to_mask_batch_idxs] // self.hop_length) - 1
                offset_idx = batch_idxs_down.min()
                max_idx = batch_idxs_up.max()

                mask = torch.arange(max_idx - offset_idx, device=device).expand(to_mask_batch_idxs.shape[0], -1)
                mask = ((batch_idxs_down - offset_idx).unsqueeze(1) <= mask) & (
                    mask < (batch_idxs_up - offset_idx).unsqueeze(1)
                )
                mask = mask.unsqueeze(-1).expand(-1, -1, self.win_length)
                masked_frames = frames[to_mask_batch_idxs, offset_idx:max_idx].masked_fill_(mask, 0)
                frames[to_mask_batch_idxs, offset_idx:max_idx] = masked_frames
        # ---

        # apply pre-emphasis first order filter on fft windows
        frames_prev = torch.roll(frames, 1, dims=-1)
        frames_prev[:, :, 0] = frames_prev[:, :, 1]
        frames = (frames - self.preemphasis * frames_prev) * 32768

        # apply fft
        S = torch.fft.rfft(fft_window * frames.view(-1, self.win_length), n=self.n_fft, dim=1)
        S = S.view(frames.shape[0], -1, S.shape[-1])
        S = S.to(torch.complex64)

        spec = torch.abs(S)
        spec_power = spec**2

        # apply triangular mel filter bank
        mel_filters = torch.from_numpy(self.mel_filters).to(device, torch.float32)
        log_spec = torch.clamp(spec_power @ mel_filters, min=1.0)
        log_spec = torch.log(log_spec)

        return log_spec