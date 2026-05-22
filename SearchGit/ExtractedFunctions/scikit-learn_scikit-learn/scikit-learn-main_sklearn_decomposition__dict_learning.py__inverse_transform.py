    def _inverse_transform(self, code, dictionary):
        """Private method allowing to accommodate both DictionaryLearning and
        SparseCoder."""
        code = check_array(code)
        # compute number of expected features in code
        expected_n_components = dictionary.shape[0]
        if self.split_sign:
            expected_n_components += expected_n_components
        if not code.shape[1] == expected_n_components:
            raise ValueError(
                "The number of components in the code is different from the "
                "number of components in the dictionary."
                f"Expected {expected_n_components}, got {code.shape[1]}."
            )
        if self.split_sign:
            n_samples, n_features = code.shape
            n_features //= 2
            code = code[:, :n_features] - code[:, n_features:]

        return code @ dictionary