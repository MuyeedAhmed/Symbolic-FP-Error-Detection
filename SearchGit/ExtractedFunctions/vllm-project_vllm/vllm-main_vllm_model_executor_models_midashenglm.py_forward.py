    def forward(self, waveform: torch.Tensor) -> torch.Tensor:
        spectrogram = F.spectrogram(
            waveform=waveform.to(torch.float32),
            pad=0,
            window=self.spectrogram_window,
            n_fft=self.config.n_fft,
            hop_length=self.config.hop_length,
            win_length=self.config.win_length,
            power=2,
            normalized=False,
            center=self.config.center,
        )
        mel_spectrogram = (spectrogram.mT @ self.melscale_fbanks.to(torch.float32)).mT
        # x has shape [batch, freq, time].
        # F.amplitude_to_DB accepts inputs shaped as:
        #   - [freq, time]
        #   - [channel, freq, time]
        #   - [..., channel, freq, time]
        # Here we insert a channel dimension of size 1 before calling it,
        # then remove that extra dimension afterward.
        log_mel_spectrogram = F.amplitude_to_DB(
            mel_spectrogram.unsqueeze(1),
            multiplier=10,
            amin=1e-10,
            db_multiplier=0,
            top_db=120,
        ).squeeze(1)
        return log_mel_spectrogram.to(waveform.dtype)