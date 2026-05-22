    def compute_whisper_melspec(
        self,
        audio_waveforms: torch.Tensor,
    ) -> torch.Tensor:
        input_dtype = audio_waveforms.dtype
        window = torch.hann_window(
            self.config.window_size, device=audio_waveforms.device
        )
        stft = torch.stft(
            audio_waveforms,
            self.config.window_size,
            self.config.hop_length,
            window=window,
            return_complex=True,
        )
        magnitudes = stft[..., :-1].abs() ** 2
        mel_spec = self.mel_filters.T @ magnitudes
        log_spec = torch.clamp(mel_spec, min=1e-10).log10()

        if global_log_mel_max := self.config.global_log_mel_max:
            if not isinstance(global_log_mel_max, float):
                raise TypeError(f"{global_log_mel_max=} needs to be of type float.")
            log_spec_max = torch.tensor(
                global_log_mel_max,
                device=log_spec.device,
                dtype=log_spec.dtype,
            )
        else:
            log_spec_max = log_spec.max()

        log_spec = torch.maximum(log_spec, log_spec_max - 8.0)
        log_spec = (log_spec + 4.0) / 4.0
        return log_spec.to(input_dtype)