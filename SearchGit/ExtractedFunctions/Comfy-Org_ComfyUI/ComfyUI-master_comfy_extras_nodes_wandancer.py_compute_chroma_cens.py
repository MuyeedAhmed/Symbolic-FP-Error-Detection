def compute_chroma_cens(y, sr=22050, hop_length=512, n_chroma=12,
                       n_octaves=7, bins_per_octave=36,
                       win_len_smooth=41, norm=2):
    """Compute Chroma Energy Normalized Statistics (CENS) features."""

    tuning = estimate_tuning(y, sr, bins_per_octave=bins_per_octave)

    fmin = 32.70319566257483  # C1 note frequency
    n_bins = n_octaves * bins_per_octave
    cqt_mag = compute_cqt(y, sr=sr, hop_length=hop_length,
                         fmin=fmin, n_bins=n_bins,
                         bins_per_octave=bins_per_octave,
                         tuning=tuning)

    chroma_map = cq_to_chroma_mapping(n_bins, bins_per_octave=bins_per_octave,
                                     n_chroma=n_chroma, fmin=fmin)
    chroma = np.dot(chroma_map, cqt_mag)

    threshold = np.finfo(chroma.dtype).tiny
    chroma_sum = np.sum(np.abs(chroma), axis=0, keepdims=True)
    chroma_sum = np.maximum(chroma_sum, threshold)
    chroma = chroma / chroma_sum

    quant_steps = [0.4, 0.2, 0.1, 0.05]
    quant_weights = [0.25, 0.25, 0.25, 0.25]
    chroma_quant = np.zeros_like(chroma)
    for step, weight in zip(quant_steps, quant_weights):
        chroma_quant += (chroma > step) * weight

    if win_len_smooth is not None and win_len_smooth > 0:
        win = scipy.signal.get_window('hann', win_len_smooth + 2, fftbins=False)
        win /= np.sum(win)
        win = win.reshape(1, -1)
        chroma_smooth = scipy.ndimage.convolve(chroma_quant, win, mode='constant')
    else:
        chroma_smooth = chroma_quant

    if norm == 2:
        threshold = np.finfo(chroma_smooth.dtype).tiny
        chroma_norm = np.sqrt(np.sum(chroma_smooth ** 2, axis=0, keepdims=True))
        chroma_norm = np.maximum(chroma_norm, threshold)
        chroma_smooth = chroma_smooth / chroma_norm
    elif norm == np.inf:
        threshold = np.finfo(chroma_smooth.dtype).tiny
        chroma_norm = np.max(np.abs(chroma_smooth), axis=0, keepdims=True)
        chroma_norm = np.maximum(chroma_norm, threshold)
        chroma_smooth = chroma_smooth / chroma_norm

    return chroma_smooth