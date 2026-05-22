def _compute_mel_spectrogram(data, sr, n_fft=2048, hop_length=512, n_mels=128):
    """Compute mel spectrogram from audio signal."""
    fft_window = scipy.signal.get_window('hann', n_fft, fftbins=True)
    if len(fft_window) < n_fft:
        lpad = int((n_fft - len(fft_window)) // 2)
        fft_window = np.pad(fft_window, (lpad, int(n_fft - len(fft_window) - lpad)), mode='constant')

    fft_window = fft_window.reshape((-1, 1))
    data_padded = np.pad(data, int(n_fft // 2), mode='constant')
    n_frames = 1 + (len(data_padded) - n_fft) // hop_length
    shape = (n_fft, n_frames)
    strides = (data_padded.strides[0], data_padded.strides[0] * hop_length)
    frames = np.lib.stride_tricks.as_strided(data_padded, shape=shape, strides=strides)

    stft_result = scipy.fft.rfft(fft_window * frames, axis=0).astype(np.complex64)
    power_spec = np.abs(stft_result) ** 2

    mel_basis = _create_mel_filterbank(sr, n_fft, n_mels=n_mels, fmin=0.0, fmax=sr / 2.0)
    mel_spec = np.dot(mel_basis, power_spec)
    return mel_spec.astype(np.float32)