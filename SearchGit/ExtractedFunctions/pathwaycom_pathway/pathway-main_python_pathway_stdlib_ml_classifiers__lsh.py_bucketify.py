    def bucketify(x: np.ndarray) -> np.ndarray:
        signs = (x @ random_hyperplanes >= 0).astype(int)
        # compute single-number bucket identifiers (i.e. single bucket int representing ANDed buckets)
        split = np.split(signs, L)
        powers = 2 ** np.arange(M).reshape(-1, 1)  # powers of two
        return np.hstack([x @ powers for x in split])