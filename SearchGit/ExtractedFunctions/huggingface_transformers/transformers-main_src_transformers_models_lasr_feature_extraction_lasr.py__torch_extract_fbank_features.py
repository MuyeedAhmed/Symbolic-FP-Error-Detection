    def _torch_extract_fbank_features(self, waveform, device="cpu"):
        # spectrogram
        window = torch.hann_window(self.win_length, periodic=False, device=device, dtype=torch.float64)
        waveform = waveform.to(torch.float64)

        # TODO: @eustlb, to be standardized
        # here we cannot use directly torch.stft because every fft frame is padded with zeros
        # due to unfold then rfft, while torch.stft unfolds with the number of fft points
        frames = waveform.unfold(-1, self.win_length, self.hop_length)
        stft = torch.fft.rfft(window * frames, n=self.n_fft)
        power_spec = torch.abs(stft) ** 2

        # log mel spectrogram
        mel_filters = self.mel_filters.to(device)
        mel_spec = torch.clamp(power_spec @ mel_filters, min=1e-5)
        mel_spec = torch.log(mel_spec)

        return mel_spec