    def make_decrypt_key(self) -> np.ndarray:
        """
        >>> hill_cipher = HillCipher(np.array([[2, 5], [1, 6]]))
        >>> hill_cipher.make_decrypt_key()
        array([[ 6, 25],
               [ 5, 26]])
        """
        det = round(np.linalg.det(self.encrypt_key))

        if det < 0:
            det = det % len(self.key_string)
        det_inv = None
        for i in range(len(self.key_string)):
            if (det * i) % len(self.key_string) == 1:
                det_inv = i
                break

        inv_key = (
            det_inv * np.linalg.det(self.encrypt_key) * np.linalg.inv(self.encrypt_key)
        )

        return self.to_int(self.modulus(inv_key))