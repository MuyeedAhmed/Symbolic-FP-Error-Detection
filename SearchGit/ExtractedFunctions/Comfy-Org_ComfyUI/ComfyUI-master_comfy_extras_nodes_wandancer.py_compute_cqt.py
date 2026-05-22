def compute_cqt(y, sr=22050, hop_length=512, fmin=None, n_bins=84, bins_per_octave=12, tuning=0.0):
    """Compute Constant-Q Transform (CQT) spectrogram."""

    def _relative_bandwidth(freqs):
        bpo = np.empty_like(freqs)
        logf = np.log2(freqs)
        bpo[0] = 1.0 / (logf[1] - logf[0])
        bpo[-1] = 1.0 / (logf[-1] - logf[-2])
        bpo[1:-1] = 2.0 / (logf[2:] - logf[:-2])
        return (2.0 ** (2.0 / bpo) - 1.0) / (2.0 ** (2.0 / bpo) + 1.0)

    def _wavelet_lengths(freqs, sr, filter_scale, alpha):
        Q = float(filter_scale) / alpha
        return Q * sr / freqs  # shape (n_bins,) floats

    def _build_wavelet(freqs_oct, sr, filter_scale, alpha_oct):
        lengths = _wavelet_lengths(freqs_oct, sr, filter_scale, alpha_oct)
        filters = []
        for ilen, freq in zip(lengths, freqs_oct):
            t = np.arange(int(-ilen // 2), int(ilen // 2), dtype=float)
            sig = (np.cos(t * 2 * np.pi * freq / sr)
                   + 1j * np.sin(t * 2 * np.pi * freq / sr)).astype(np.complex64)
            sig *= scipy.signal.get_window('hann', len(sig), fftbins=True)
            l1 = np.sum(np.abs(sig))
            tiny = np.finfo(np.float32).tiny
            sig /= max(l1, tiny)
            filters.append(sig)
        max_len = max(lengths)
        n_fft = int(2.0 ** np.ceil(np.log2(max_len)))
        out = np.zeros((len(filters), n_fft), dtype=np.complex64)
        for k, f in enumerate(filters):
            lpad = int((n_fft - len(f)) // 2)
            out[k, lpad: lpad + len(f)] = f
        return out, lengths

    def _resample_half(y):
        ratio = 0.5
        n_samples = int(np.ceil(len(y) * ratio))
        # Kaiser-windowed FIR matches librosa/soxr more closely than scipy's default Hamming filter
        L = 2
        h = scipy.signal.firwin(160 * L + 1, 0.96 / L, window=('kaiser', 6.5))
        y_hat = scipy.signal.resample_poly(y.astype(np.float32), 1, 2, window=h)
        if len(y_hat) > n_samples:
            y_hat = y_hat[:n_samples]
        elif len(y_hat) < n_samples:
            y_hat = np.pad(y_hat, (0, n_samples - len(y_hat)))
        y_hat /= np.sqrt(ratio)
        return y_hat.astype(np.float32)

    def _sparsify_rows(x, quantile=0.01):
        mags = np.abs(x)
        norms = np.sum(mags, axis=1, keepdims=True)
        norms = np.where(norms == 0, 1.0, norms)
        mag_sort = np.sort(mags, axis=1)
        cumulative_mag = np.cumsum(mag_sort / norms, axis=1)
        threshold_idx = np.argmin(cumulative_mag < quantile, axis=1)
        x_sparse = scipy.sparse.lil_matrix(x.shape, dtype=x.dtype)
        for i, j in enumerate(threshold_idx):
            idx = np.where(mags[i] >= mag_sort[i, j])
            x_sparse[i, idx] = x[i, idx]
        return x_sparse.tocsr()

    if fmin is None:
        fmin = 32.70319566257483  # C1 note frequency

    fmin = fmin * (2.0 ** (tuning / bins_per_octave))
    freqs = fmin * (2.0 ** (np.arange(n_bins) / bins_per_octave))

    alpha = _relative_bandwidth(freqs)
    lengths = _wavelet_lengths(freqs, float(sr), 1, alpha)

    n_octaves = int(np.ceil(float(n_bins) / bins_per_octave))
    n_filters = min(bins_per_octave, n_bins)

    cqt_resp = []
    my_y = y.astype(np.float32)
    my_sr = float(sr)
    my_hop = int(hop_length)

    for i in range(n_octaves):
        if i == 0:
            sl = slice(-n_filters, None)
        else:
            sl = slice(-n_filters * (i + 1), -n_filters * i)

        freqs_oct = freqs[sl]
        alpha_oct = alpha[sl]

        basis, basis_lengths = _build_wavelet(freqs_oct, my_sr, 1, alpha_oct)
        n_fft_oct = basis.shape[1]

        # Frequency-domain normalisation
        basis = basis.astype(np.complex64)
        basis *= basis_lengths[:, np.newaxis] / float(n_fft_oct)
        fft_basis = scipy.fft.fft(basis, n=n_fft_oct, axis=1)[:, :(n_fft_oct // 2) + 1]
        fft_basis = _sparsify_rows(fft_basis, quantile=0.01)
        fft_basis = fft_basis * np.sqrt(sr / my_sr)

        y_pad = np.pad(my_y, int(n_fft_oct // 2), mode='constant')
        n_frames = 1 + (len(y_pad) - n_fft_oct) // my_hop
        frames = np.lib.stride_tricks.as_strided(
            y_pad,
            shape=(n_fft_oct, n_frames),
            strides=(y_pad.strides[0], y_pad.strides[0] * my_hop),
        )
        stft_result = scipy.fft.rfft(frames, axis=0)
        cqt_resp.append(fft_basis.dot(stft_result))

        if my_hop % 2 == 0:
            my_hop //= 2
            my_sr /= 2.0
            my_y = _resample_half(my_y)

    max_col = min(c.shape[-1] for c in cqt_resp)
    cqt_out = np.empty((n_bins, max_col), dtype=np.complex64)
    end = n_bins
    for c_i in cqt_resp:
        n_oct = c_i.shape[0]
        if end < n_oct:
            cqt_out[:end, :] = c_i[-end:, :max_col]
        else:
            cqt_out[end - n_oct:end, :] = c_i[:, :max_col]
        end -= n_oct

    cqt_out /= np.sqrt(lengths)[:, np.newaxis]
    return np.abs(cqt_out).astype(np.float32)