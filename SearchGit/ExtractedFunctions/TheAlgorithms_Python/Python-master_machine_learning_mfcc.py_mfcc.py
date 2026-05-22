def mfcc(
    audio: np.ndarray,
    sample_rate: int,
    ftt_size: int = 1024,
    hop_length: int = 20,
    mel_filter_num: int = 10,
    dct_filter_num: int = 40,
) -> np.ndarray:
    """
    Calculate Mel Frequency Cepstral Coefficients (MFCCs) from an audio signal.

    Args:
        audio: The input audio signal.
        sample_rate: The sample rate of the audio signal (in Hz).
        ftt_size: The size of the FFT window (default is 1024).
        hop_length: The hop length for frame creation (default is 20ms).
        mel_filter_num: The number of Mel filters (default is 10).
        dct_filter_num: The number of DCT filters (default is 40).

    Returns:
        A matrix of MFCCs for the input audio.

    Raises:
        ValueError: If the input audio is empty.

    Example:
    >>> sample_rate = 44100  # Sample rate of 44.1 kHz
    >>> duration = 2.0  # Duration of 1 second
    >>> t = np.linspace(0, duration, int(sample_rate * duration), endpoint=False)
    >>> audio = 0.5 * np.sin(2 * np.pi * 440.0 * t)  # Generate a 440 Hz sine wave
    >>> mfccs = mfcc(audio, sample_rate)
    >>> mfccs.shape
    (40, 101)
    """
    logging.info(f"Sample rate: {sample_rate}Hz")
    logging.info(f"Audio duration: {len(audio) / sample_rate}s")
    logging.info(f"Audio min: {np.min(audio)}")
    logging.info(f"Audio max: {np.max(audio)}")

    # normalize audio
    audio_normalized = normalize(audio)

    logging.info(f"Normalized audio min: {np.min(audio_normalized)}")
    logging.info(f"Normalized audio max: {np.max(audio_normalized)}")

    # frame audio into
    audio_framed = audio_frames(
        audio_normalized, sample_rate, ftt_size=ftt_size, hop_length=hop_length
    )

    logging.info(f"Framed audio shape: {audio_framed.shape}")
    logging.info(f"First frame: {audio_framed[0]}")

    # convert to frequency domain
    # For simplicity we will choose the Hanning window.
    window = get_window("hann", ftt_size, fftbins=True)
    audio_windowed = audio_framed * window

    logging.info(f"Windowed audio shape: {audio_windowed.shape}")
    logging.info(f"First frame: {audio_windowed[0]}")

    audio_fft = calculate_fft(audio_windowed, ftt_size)
    logging.info(f"fft audio shape: {audio_fft.shape}")
    logging.info(f"First frame: {audio_fft[0]}")

    audio_power = calculate_signal_power(audio_fft)
    logging.info(f"power audio shape: {audio_power.shape}")
    logging.info(f"First frame: {audio_power[0]}")

    filters = mel_spaced_filterbank(sample_rate, mel_filter_num, ftt_size)
    logging.info(f"filters shape: {filters.shape}")

    audio_filtered = np.dot(filters, np.transpose(audio_power))
    audio_log = 10.0 * np.log10(audio_filtered)
    logging.info(f"audio_log shape: {audio_log.shape}")

    dct_filters = discrete_cosine_transform(dct_filter_num, mel_filter_num)
    cepstral_coefficents = np.dot(dct_filters, audio_log)

    logging.info(f"cepstral_coefficents shape: {cepstral_coefficents.shape}")
    return cepstral_coefficents