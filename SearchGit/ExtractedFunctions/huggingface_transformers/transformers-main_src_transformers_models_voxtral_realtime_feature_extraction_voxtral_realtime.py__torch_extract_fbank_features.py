    def _torch_extract_fbank_features(self, waveform, device: str = "cpu", center: bool = True):
        window = torch.hann_window(self.n_fft, device=device)
        stft = torch.stft(waveform, self.n_fft, self.hop_length, window=window, return_complex=True, center=center)
        magnitudes = stft[..., :-1].abs() ** 2

        mel_filters = torch.from_numpy(self.mel_filters).to(device, torch.float32)
        mel_spec = mel_filters.T @ magnitudes

        log_spec = torch.clamp(mel_spec, min=1e-10).log10()
        if self.global_log_mel_max is not None:
            log_spec_max = torch.tensor(
                self.global_log_mel_max,
                device=log_spec.device,
                dtype=log_spec.dtype,
            )
        else:
            log_spec_max = log_spec.max()

        log_spec = torch.maximum(log_spec, log_spec_max - 8.0)
        log_spec = (log_spec + 4.0) / 4.0
        if device != "cpu":
            log_spec = log_spec.detach().cpu()
        return log_spec