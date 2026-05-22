    def _apply_mel_filters(
        self, stft_output: torch.Tensor, mel_filters: torch.Tensor
    ) -> torch.Tensor:
        magnitudes = stft_output.real.square() + stft_output.imag.square()
        mel_spec = mel_filters @ magnitudes
        mel_spec = torch.log(mel_spec + LOG_ZERO_GUARD_VALUE)
        return mel_spec.permute(0, 2, 1)